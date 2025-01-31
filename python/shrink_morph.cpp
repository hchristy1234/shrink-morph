#include "functions.h"
#include "newton.h"
#include "parameterization.h"
#include "path_extraction.h"
#include "simulation_utils.h"
#include "stretch_angles.h"
#include "stripe_patterns.h"
#include "timer.h"
#include "SGNSolver.h"

#include <Eigen/Dense>
#include <igl/loop.h>
#include <igl/readOBJ.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
// #include "geometrycentral/surface/stripe_patterns.h"

#include <string>
#include <tuple>

namespace nb = nanobind;
using namespace nb::literals;

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi> readFromOBJ(std::string fileName)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(!igl::readOBJ(fileName, V, F))
  {
    std::cout << "File " << fileName << " not found\n";
  }
  V.col(1) *= -1;
  return std::make_tuple(V, F);
}

// std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXi>
// subdivide(const nb::DRef<Eigen::MatrixXd>& V, const nb::DRef<Eigen::MatrixXd>& P, const nb::DRef<Eigen::MatrixXi>& F)
// {
//   Eigen::VectorXd x;
//   const auto& [NV, NP, NF, nx] = subdivide(V, P, F, x);

//   return std::make_tuple(NV, NP, NF);
// }

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXd>
subdivide(const nb::DRef<Eigen::MatrixXd>& V,
          const nb::DRef<Eigen::MatrixXd>& P,
          const nb::DRef<Eigen::MatrixXi>& F,
          const nb::DRef<Eigen::VectorXd>& x,
          double targetLength)
{
  Eigen::MatrixXi NF;
  Eigen::SparseMatrix<double> S;

  igl::loop(V.rows(), F, S, NF);

  Eigen::MatrixXd NV = S * V;
  Eigen::MatrixXd NP = S * P;
  Eigen::MatrixXd nx;
  if(x.size() == S.cols())
    nx = S * x;

  if(targetLength > 0)
  {
    // create geometry-central objects
    using namespace geometrycentral::surface;
    ManifoldSurfaceMesh mesh(NF);
    VertexPositionGeometry geometry(mesh, NV);

    double avgEdgeLength = 0;
    geometry.requireVertexPositions();
    for(Edge e: mesh.edges())
      avgEdgeLength += norm(geometry.vertexPositions[e.firstVertex()] - geometry.vertexPositions[e.secondVertex()]);
    avgEdgeLength /= mesh.nEdges();

    while(avgEdgeLength > targetLength)
    {
      Eigen::MatrixXi tempF = NF;
      igl::loop(NV.rows(), tempF, S, NF);
      NV = S * NV;
      NP = S * NP;
      if(nx.size() == S.cols())
        nx = S * nx;

      avgEdgeLength /= 2;
    }
  }
  return std::make_tuple(NV, NP, NF, nx);
}

void simulation(nb::DRef<Eigen::MatrixXd> V,
                nb::DRef<Eigen::MatrixXd> P,
                const nb::DRef<Eigen::MatrixXi>& F,
                const nb::DRef<Eigen::VectorXd>& theta2,
                double E1,
                double lambda1,
                double lambda2,
                // double deltaLambda,
                double thickness,
                double width,
                int n_iter,
                double lim)
{
  Timer timer("Simulation");

  using namespace geometrycentral::surface;

  const double scaleFactor = width / (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
  V *= scaleFactor;
  P *= scaleFactor;

  // create geometry-central objects
  ManifoldSurfaceMesh mesh(F);

  // Find face closest to mesh center and fix its vertices
  std::vector<int> fixedIdx = findCenterFaceIndices(P, F);

  FaceData<Eigen::Matrix2d> MrInv = precomputeSimData(mesh, P, F);
  FaceData<double> theta1 = computeStretchAngles(mesh, V, F, MrInv);
  VertexData<double> _theta2(mesh, theta2);

  // Define simulation function
  auto func = simulationFunction(mesh, MrInv, theta1, _theta2, E1, lambda1, lambda2, 0, thickness);

  // ---------------- (Projected) Newton optimization ----------------
  LLTSolver solver;
  Eigen::VectorXd x = V.reshaped<Eigen::RowMajor>();

  // Newton algorithm
  newton(x, func, solver, n_iter, lim, true, fixedIdx);

  V = x.reshaped<Eigen::RowMajor>(V.rows(), 3);
}

Eigen::VectorXd directionsOptimization(nb::DRef<Eigen::MatrixXd> V,
                                       nb::DRef<Eigen::MatrixXd> targetV,
                                       nb::DRef<Eigen::MatrixXd> P,
                                       const Eigen::MatrixXi& F,
                                       double E1,
                                       double lambda1,
                                       double lambda2,
                                      //  double deltaLambda,
                                       double thickness,
                                       double width,
                                       int n_iter,
                                       double lim,
                                       double wM,
                                       double wL)
{

  Timer timer("Directions optimization");

  using namespace geometrycentral::surface;

  const double scaleFactor = width / (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
  V *= scaleFactor;
  targetV *= scaleFactor;
  P *= scaleFactor;

  // create geometry-central objects
  ManifoldSurfaceMesh mesh(F);
  VertexPositionGeometry geometry(mesh, V);

  // Find face closest to mesh center and fix its vertices
  std::vector<int> fixedIdx = findCenterFaceIndices(P, F);

  FaceData<Eigen::Matrix2d> MrInv = precomputeSimData(mesh, P, F);
  FaceData<double> theta1 = computeStretchAngles(mesh, V, F, MrInv);
  VertexData<double> theta2(mesh, 0);

  // Define optimization function
  auto adjointFunc = adjointFunction(geometry, F, MrInv, theta1, E1, lambda1, lambda2, 0, thickness);

  // Optimize this energy function using SGN [Zehnder et al. 2021]
  V = sparse_gauss_newton(geometry, targetV, MrInv, theta1, theta2, adjointFunc, fixedIdx, n_iter, lim, wM, wL, E1,
                          lambda1, lambda2, 0, thickness);

  return theta2.toVector();
}

// std::vector<Eigen::MatrixXd> generateTrajectories(const nb::DRef<Eigen::MatrixXd>& P,
//                                                   const nb::DRef<Eigen::MatrixXi>& F,
//                                                   const nb::DRef<Eigen::MatrixXd>& directions,
//                                                   double spacing)
// {
//   using namespace geometrycentral;
//   using namespace geometrycentral::surface;

//   // create geometry-central objects
//   ManifoldSurfaceMesh mesh(F);
//   Eigen::MatrixXd P_3D(P.rows(), 3);
//   P_3D.leftCols(2) = P;
//   P_3D.col(2).setZero();
//   VertexPositionGeometry geometryUV(mesh, P_3D);

//   geometryUV.requireVertexTangentBasis();
//   VertexData<double> frequencies(mesh, 1.0 / spacing);

//   // compute direction field
//   VertexData<Vector2> directionField(mesh);
//   for(size_t i = 0; i < mesh.nVertices(); ++i)
//   {
//     // interpolate orientations w.r.t layer
//     directionField[i] = {directions(i, 0), directions(i, 1)};
//     // express vectors in their respective vertex local bases (the stripes algorithm expects that)
//     auto basisVector = geometryUV.vertexTangentBasis[i][0];
//     Vector2 base{basisVector.x, basisVector.y};
//     directionField[i] = -directionField[i].pow(2) / base.pow(2);
//   }

//   const auto& [stripeValues, stripeIndices, fieldIndices] = computeStripePattern(geometryUV, frequencies, directionField);
//   const auto& [vertices, edges] = extractPolylinesFromStripePattern(geometryUV, stripeValues, stripeIndices, fieldIndices, directionField, true);

//   auto polylines = edgeToPolyline(vertices, edges);
//   polylines = orderPolylines(polylines);
//   return simplifyPolylines(polylines);
// }

struct StripeAlgo
{
  Eigen::SparseMatrix<double> massMatrix;
  Eigen::VectorXd u;
  LDLTSolver solver;
  // geometrycentral::surface::VertexPositionGeometry geometryUV;

  StripeAlgo(const nb::DRef<Eigen::MatrixXd>& P, const nb::DRef<Eigen::MatrixXi>& F)
  {
    using namespace geometrycentral;
    using namespace geometrycentral::surface;

    ManifoldSurfaceMesh mesh(F);
    Eigen::MatrixXd P_3D(P.rows(), 3);
    P_3D.leftCols(2) = P;
    P_3D.col(2).setZero();
    VertexPositionGeometry geometryUV(mesh, P_3D);

    massMatrix = computeRealVertexMassMatrix(geometryUV);
    u = Eigen::VectorXd::Random(mesh.nVertices() * 2);
  }

  std::vector<Eigen::MatrixXd> generateTrajectory(const nb::DRef<Eigen::MatrixXd>& P,
                                                  const nb::DRef<Eigen::MatrixXi>& F,
                                                  const nb::DRef<Eigen::VectorXd>& theta1,
                                                  const nb::DRef<Eigen::VectorXd>& theta2,
                                                  double layerHeight,
                                                  double spacing,
                                                  int nLayers,
                                                  int i)
  {
    using namespace geometrycentral;
    using namespace geometrycentral::surface;

    // create geometry-central objects
    ManifoldSurfaceMesh mesh(F);
    Eigen::MatrixXd P_3D(P.rows(), 3);
    P_3D.leftCols(2) = P;
    P_3D.col(2).setZero();
    VertexPositionGeometry geometryUV(mesh, P_3D);

    return generateOneLayer(geometryUV, theta1, theta2, massMatrix, u, solver, i, nLayers, layerHeight, spacing);
  }
};

Eigen::VectorXd vertexBasedStretchAngles(const nb::DRef<Eigen::MatrixXd>& V,
                                         const nb::DRef<Eigen::MatrixXd>& P,
                                         const nb::DRef<Eigen::MatrixXi>& F)
{
  using namespace geometrycentral;
  using namespace geometrycentral::surface;

  // create geometry-central objects
  ManifoldSurfaceMesh mesh(F);
  FaceData<double> theta1 = computeStretchAngles(mesh, V, P, F);

  // convert face angles into vertex angles
  Eigen::VectorXd vertexAngles(V.rows());
  for(int i = 0; i < V.rows(); ++i)
  {
    Vector2 direction = {0, 0};
    Vertex v = mesh.vertex(i);
    for(Face f: v.adjacentFaces())
    {
      if(!f.isBoundaryLoop())
      {
        direction += Vector2::fromAngle(2 * theta1[f]);
      }
    }
    direction = direction.normalize().pow(0.5);
    vertexAngles(i) = direction.arg();
  }

  return vertexAngles;
}

NB_MODULE(shrink_morph_py, m)
{
  m.def("read_from_OBJ", &readFromOBJ);
  m.def("subdivide", &subdivide, "V"_a, "P"_a, "F"_a, "x"_a, "target_edge_length"_a = 0);
  m.def("compute_SVD_data", [](const nb::DRef<Eigen::MatrixXd>& V, const nb::DRef<Eigen::MatrixXd>& P,
                               const nb::DRef<Eigen::MatrixXi>& F) { return computeSVDdata(V, P, F); });
  m.def("parameterization",
        [](const nb::DRef<Eigen::MatrixXd>& V, Eigen::MatrixXi& F, double lambda1, double lambda2, double wD,
           int n_iter, double lim) { return parameterization(V, F, lambda1, lambda2, wD, n_iter, lim); });
  m.def("reparameterization",
        [](const nb::DRef<Eigen::MatrixXd>& V, Eigen::MatrixXd& P, const nb::DRef<Eigen::MatrixXi>& F, double lambda1,
           double lambda2, double wD, int n_iter, double lim) {
          parameterization(V, P, F, lambda1, lambda2, wD, n_iter, lim);
          return P;
        });
  m.def("simulation", &simulation);
  m.def("directions_optimization", &directionsOptimization);
  // m.def("generate_trajectories", &generateTrajectories);
  m.def("vertex_based_stretch_angles", &vertexBasedStretchAngles);
  nb::class_<StripeAlgo>(m, "StripeAlgo")
      .def(nb::init<const nb::DRef<Eigen::MatrixXd>&, const nb::DRef<Eigen::MatrixXi>&>())
      .def("generate_one_layer", &StripeAlgo::generateTrajectory);

  nb::class_<SGNSolver>(m, "SGNSolver")
      .def(nb::init<const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                    double, double, double, double>())
      .def("solve_one_step", &SGNSolver::solveOneStep)
      .def("optimizedV", &SGNSolver::vertices)
      .def("distance", &SGNSolver::distance)
      .def("decrement", &SGNSolver::newton_decrement);
}
