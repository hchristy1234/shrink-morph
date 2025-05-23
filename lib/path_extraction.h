#pragma once

#include "solvers.h"

#include <geometrycentral/surface/vertex_position_geometry.h>

std::vector<std::vector<geometrycentral::Vector3>>
orderPolylines(const std::vector<std::vector<geometrycentral::Vector3>>& isolines);

std::vector<Eigen::MatrixXd> simplifyPolylines(const std::vector<std::vector<geometrycentral::Vector3>>& polylines);

std::vector<std::vector<Eigen::MatrixXd>> generatePaths(geometrycentral::surface::EmbeddedGeometryInterface& geometry,
                                                        const Eigen::VectorXd& theta1,
                                                        const Eigen::VectorXd& theta2,
                                                        int nLayers,
                                                        double spacing);

std::vector<Eigen::MatrixXd> generateOneLayer(geometrycentral::surface::EmbeddedGeometryInterface& geometry,
                                              const Eigen::VectorXd& theta,
                                              const Eigen::SparseMatrix<double>& massMatrix,
                                              Eigen::VectorXd& u,
                                              LDLTSolver& solver,
                                              bool patternAnalyzed,
                                              double spacing);

void writePaths(const std::string& filename, const std::vector<Eigen::MatrixXd>& paths, double height);

std::vector<std::vector<geometrycentral::Vector3>> edgeToPolyline(const std::vector<geometrycentral::Vector3>& points,
                                                                  const std::vector<std::array<size_t, 2>>& edges);

std::tuple<geometrycentral::surface::VertexData<double>,
           geometrycentral::surface::FaceData<double>,
           geometrycentral::surface::EdgeData<double>>
hodgeDecomposition(geometrycentral::surface::VertexPositionGeometry& geometry,
                   const geometrycentral::surface::EdgeData<double>& oneForm);

Eigen::VectorXd
curlFreeParameterization(geometrycentral::surface::VertexPositionGeometry& geometry,
                         const geometrycentral::surface::FaceData<geometrycentral::Vector3>& vectorField);

Eigen::VectorXd
curlFreeParameterization(geometrycentral::surface::VertexPositionGeometry& geometry,
                         const geometrycentral::surface::FaceData<geometrycentral::Vector3>& vectorField,
                         LDLTSolver& solver);

Eigen::VectorXd
curlFreeParameterization(const Eigen::MatrixXd& P, const Eigen::MatrixXi& F, const Eigen::VectorXd& theta);

std::vector<std::vector<geometrycentral::Vector3>>
orderPolylinesNew(const std::vector<std::vector<geometrycentral::Vector3>>& isolines,
                  const Eigen::MatrixX2d& P,
                  const Eigen::MatrixXi& F,
                  const Eigen::VectorXd& theta);

std::vector<std::vector<geometrycentral::Vector3>>
orderPolylinesNew(const std::vector<std::vector<geometrycentral::Vector3>>& isolines,
                  const Eigen::MatrixX2d& P,
                  const Eigen::VectorXd& vCoordinate);