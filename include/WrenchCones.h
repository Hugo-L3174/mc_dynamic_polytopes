#pragma once
#include <mc_control/fsm/Controller.h>
// #include <eigen-cdd/Polyhedron.h>

Eigen::MatrixXd linearizedFrictionCone(int numberOfFrictionSides, Eigen::Matrix3d m_rotation, double m_frictionCoef);

// Compute the directions for the generators of a linearized friction cone from its orientation in world frame and the
// precision of the linearization
std::vector<Eigen::Vector3d> generatePolyhedralConeGens(int numberOfFrictionSides,
                                                        Eigen::Matrix3d m_rotation,
                                                        double m_frictionCoef,
                                                        double scale);

// Compute the 6 x n sized matrix of the generators for a single contact
// application point is around what the generators must be applied
Eigen::MatrixXd compute6DGeneratorsMatrixSingleCone(
    Eigen::Vector3d applicationPoint,
    int numberOfFrictionSides,
    std::pair<std::pair<double, double>, sva::PTransformd> contactSurface,
    double m_frictionCoef);

// Compute the concatenated 6 x n sized matrix of all the generators for the given vector of contacts (redundant
// description but if R-rep, ie unbounded cones then this is already the R-rep of the mink sum)
Eigen::MatrixXd compute6DGeneratorsMatrixRaysCones(
    Eigen::Vector3d applicationPoint,
    int numberOfFrictionSides,
    std::vector<std::pair<std::pair<double, double>, sva::PTransformd>> contactSurfaces,
    double m_frictionCoef);

// Compute the 3 x n sized matrix of the generators for a single contact
// application point is around what the generators must be applied
Eigen::MatrixXd compute3DGeneratorsMatrixSingleCone(
    Eigen::Vector3d applicationPoint,
    int numberOfFrictionSides,
    std::pair<std::pair<double, double>, sva::PTransformd> contactSurface,
    double m_frictionCoef);

// Compute the concatenated 3 x n sized matrix of all the generators for the given vector of contacts (redundant
// description but if R-rep, ie unbounded cones then this is already the R-rep of the mink sum)
Eigen::MatrixXd compute3DGeneratorsMatrixRaysCones(
    Eigen::Vector3d applicationPoint,
    int numberOfFrictionSides,
    std::vector<std::pair<std::pair<double, double>, sva::PTransformd>> contactSurfaces,
    double m_frictionCoef);

Eigen::Matrix3d skewMatrix(const Eigen::Vector3d v);
