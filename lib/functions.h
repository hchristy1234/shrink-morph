#pragma once

#include <TinyAD/Support/GeometryCentral.hh>
#include <TinyAD/ScalarFunction.hh>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>

TinyAD::ScalarFunction<3, double, geometrycentral::surface::VertexRangeF::Etype>
simulationFunction(geometrycentral::surface::SurfaceMesh& mesh,
                   const geometrycentral::surface::FaceData<Eigen::Matrix2d>& MrInv,
                   const geometrycentral::surface::FaceData<double>& theta1,
                   const geometrycentral::surface::VertexData<double>& theta2,
                   double E1,
                   double lambda1,
                   double lambda2,
                   double deltaLambda,
                   double thickness,
                   double E2 = 1);

TinyAD::ScalarFunction<1, double, Eigen::Index>
adjointFunction(geometrycentral::surface::IntrinsicGeometryInterface& geometry,
                const Eigen::MatrixXi& F,
                const geometrycentral::surface::FaceData<Eigen::Matrix2d>& MrInv,
                const geometrycentral::surface::FaceData<double>& theta1,
                double E1,
                double lambda1,
                double lambda2,
                double deltaLambda,
                double thickness,
                double E2 = 1);
