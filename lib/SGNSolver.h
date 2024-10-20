#include "solvers.h"

#include <Eigen/Core>
#include <TinyAD/Support/GeometryCentral.hh>
#include <TinyAD/ScalarFunction.hh>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

class SGNSolver
{
public:
  SGNSolver(const Eigen::MatrixXd& targetV,
            const Eigen::MatrixXd& P,
            const Eigen::MatrixXi& F,
            double _E1,
            double _lambda1,
            double _lambda2,
            double _deltaLambda,
            double _thickness);

  double distance(const Eigen::VectorXd& _theta);
  Eigen::VectorXd distanceGrad(const Eigen::VectorXd& _theta);
  std::tuple<Eigen::VectorXd, Eigen::VectorXd> solveOneStep();

  double newton_decrement() { return decrement; }
  Eigen::MatrixXd vertices() { return x.reshaped<Eigen::RowMajor>(); }
  Eigen::VectorXd thetas() { return theta; }

private:
  LUSolver solver;
  LLTSolver adjointSolver;
  geometrycentral::surface::ManifoldSurfaceMesh mesh;
  geometrycentral::surface::VertexPositionGeometry geometry;
  geometrycentral::surface::FaceData<double> theta1;
  geometrycentral::surface::FaceData<Eigen::Matrix2d> MrInv;

  TinyAD::ScalarFunction<1, double, Eigen::Index> adjointFunc;

  Eigen::SparseMatrix<double> P;
  Eigen::SparseMatrix<double> H;
  Eigen::SparseMatrix<double> HGN;
  Eigen::SparseMatrix<double> M_theta;
  Eigen::SparseMatrix<double> L;

  Eigen::VectorXd x;
  Eigen::VectorXd theta;
  Eigen::VectorXd xTarget;
  Eigen::VectorXd masses;
  Eigen::MatrixXi F;

  std::vector<int> fixedIdx;
  int max_iters;
  double lim;
  double wM;
  double wL;
  double E1;
  double lambda1;
  double lambda2;
  double deltaLambda;
  double thickness;
  double decrement;
  int iter;
};