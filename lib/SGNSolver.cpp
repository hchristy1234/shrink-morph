#include "SGNSolver.h"

#include "functions.h"
#include "newton.h"
#include "timer.h"
#include "parameterization.h"
#include "simulation_utils.h"
#include "stretch_angles.h"

#include <TinyAD/Utils/NewtonDecrement.hh>

using namespace geometrycentral::surface;

SGNSolver::SGNSolver(const Eigen::MatrixXd& targetV,
                     const Eigen::MatrixXd& _P,
                     const Eigen::MatrixXi& _F,
                     double _E1,
                     double _lambda1,
                     double _lambda2,
                     double _deltaLambda,
                     double _thickness)
    : wM(0.01),
      wL(0.01),
      lim(1e-6),
      max_iters(1000),
      E1(_E1),
      lambda1(_lambda1),
      lambda2(_lambda2),
      deltaLambda(_deltaLambda),
      F(_F),
      mesh(_F),
      geometry(mesh, targetV),
      thickness(_thickness),
      iter(0)
{
  MrInv = precomputeSimData(mesh, _P, F);
  theta1 = computeStretchAngles(mesh, targetV, F, MrInv);
  
  adjointFunc = adjointFunction(geometry, F, MrInv, theta1, E1, lambda1, lambda2, deltaLambda, thickness);


  fixedIdx = findCenterFaceIndices(_P, F);

  geometry.requireCotanLaplacian();
  geometry.requireFaceAreas();
  geometry.requireVertexIndices();

  // build vector of Voronoi areas for the distance metric
  masses = Eigen::VectorXd::Zero(targetV.size());

  // divide M by the total mesh area
  double totalArea = 0;
  for(Face f: mesh.faces())
  {
    for(Vertex v: f.adjacentVertices())
    {
      masses(3 * geometry.vertexIndices[v]) += geometry.faceAreas[f] / 3.;
      masses(3 * geometry.vertexIndices[v] + 1) += geometry.faceAreas[f] / 3.;
      masses(3 * geometry.vertexIndices[v] + 2) += geometry.faceAreas[f] / 3.;
    }
    totalArea += geometry.faceAreas[f];
  }
  masses /= totalArea;

  L = geometry.cotanLaplacian;

  // Mass matrix theta
  M_theta = Eigen::SparseMatrix<double>(targetV.rows(), targetV.rows());
  M_theta.reserve(targetV.rows());
  for(int i = 0; i < targetV.rows(); ++i)
    M_theta.insert(i, i) = totalArea * masses(3 * i);

  theta = Eigen::VectorXd::Zero(targetV.rows());
  xTarget = targetV.reshaped<Eigen::RowMajor>();
  x = xTarget;

  // Build matrix P
  P = projectionMatrix(fixedIdx, x.size());

  // Hessian matrix H
  Eigen::VectorXd X(targetV.size() + theta.size());
  X.head(targetV.size()) = x;
  X.tail(theta.size()) = theta;
  H = adjointFunc.eval_hessian(X);

  // Build HGN matrix
  HGN = buildHGN(2 * masses, P, 2 * wM * M_theta + 2 * wL * L, H);
}

double SGNSolver::distance(const Eigen::VectorXd& th)
{
  VertexData<double> theta2(mesh, th);
  auto simFunc = simulationFunction(mesh, MrInv, theta1, theta2, E1, lambda1, lambda2, deltaLambda, thickness);
  newton(x, simFunc, adjointSolver, max_iters, lim, false, fixedIdx);

  return (x - xTarget).dot(masses.cwiseProduct(x - xTarget)) + wM * th.dot(M_theta * th) + wL * th.dot(L * th);
}

Eigen::VectorXd SGNSolver::distanceGrad(const Eigen::VectorXd& th)
{
  Eigen::VectorXd X(x.size() + th.size());
  X.head(x.size()) = x;
  X.tail(th.size()) = th;
  H = adjointFunc.eval_hessian(X);

  for(int j = 0; j < x.size(); ++j)
    H.coeffRef(j, j) += 1e-10;

  Eigen::SparseMatrix<double> A = (P * H.block(0, 0, x.size(), x.size()) * P.transpose()).eval();

  adjointSolver.factorize(A);
  if(adjointSolver.info() != Eigen::Success)
  {
    auto [f, g, A_proj] = adjointFunc.eval_with_hessian_proj(X);
    A_proj = (P * A_proj.block(0, 0, x.size(), x.size()) * P.transpose()).eval();

    A = 0.9 * A + 0.1 * A_proj;
    adjointSolver.factorize(A);
    if(adjointSolver.info() != Eigen::Success)
      adjointSolver.factorize(A_proj);
  }

  Eigen::VectorXd b = P * masses.cwiseProduct(x - xTarget);
  Eigen::VectorXd dir = adjointSolver.solve(b);
  if(adjointSolver.info() != Eigen::Success)
    std::cout << "Solver error\n";

  dir = P.transpose() * dir;

  return -2 * H.block(x.size(), 0, th.size(), x.size()) * dir + 2 * wM * M_theta * th + 2 * wL * L * th;
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd> SGNSolver::solveOneStep()
{
  double f = distance(theta);
  Eigen::VectorXd g = distanceGrad(theta);

  Eigen::VectorXd b(2 * x.size() - 2 * fixedIdx.size() + theta.size());
  b.setZero();
  b.segment(x.size() - fixedIdx.size(), theta.size()) = -g;

  // Update HGN
  updateHGN(HGN, P, H);

  if(iter == 0)
    solver.compute(HGN);
  else
    solver.factorize(HGN);

  if(solver.info() != Eigen::Success)
  {
    std::cout << "Solver error\n";
    return std::tie(x, theta);
  }

  Eigen::VectorXd d = solver.solve(b);
  Eigen::VectorXd deltaTheta = d.segment(x.size() - fixedIdx.size(), theta.size());
  Eigen::VectorXd deltaX = d.segment(0, x.size() - fixedIdx.size());
  deltaX = P.transpose() * deltaX;

  // LINE SEARCH
  Eigen::VectorXd x_old = x;
  double s = lineSearch(
      theta, deltaTheta, f, g, [this](const Eigen::VectorXd& th) { return distance(th); },
      [&](double s) { x = x_old + s * deltaX; });
  if(s < 0)
  {
    std::cout << "Line search failed\n";
    return std::tie(x, theta);
  }
  theta += s * deltaTheta;

  decrement = TinyAD::newton_decrement(deltaTheta, g);

  std::cout << "Decrement in iteration " << iter << ": " << decrement
            << "\tDistance: " << (x - xTarget).dot(masses.cwiseProduct(x - xTarget)) << "\tStep size: " << s
            << std::endl;
  // if(TinyAD::newton_decrement(deltaTheta, g) < lim || solver.info() != Eigen::Success)
  //   break;

  iter++;

  return std::tie(x, theta);
}