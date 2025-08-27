#pragma once

#include <mc_rtc/gui.h>
#include <politopix/politopixAPI.h>

namespace mc_dynamic_polytopes
{

// Update faces vector for polytope display with origin at contact pose (no offset by default)
void update3DPolyTrianglesPolitopix(boost::shared_ptr<Polytope_Rn> & polytope,
                                    std::vector<std::array<Eigen::Vector3d, 3>> & resultTriangles,
                                    double guiScale,
                                    sva::PTransformd contactPose = sva::PTransformd(Eigen::Vector3d{0.0, 0.0, 0.0}));

// Updates the faces vector of the 3D force and 3D moment in the case of a 6D polytope
void update6DPolyTrianglesPolitopix(boost::shared_ptr<Polytope_Rn> & polytope,
                                    std::vector<std::array<Eigen::Vector3d, 3>> & resultMomentTriangles,
                                    std::vector<std::array<Eigen::Vector3d, 3>> & resultForceTriangles,
                                    double guiScale,
                                    sva::PTransformd contactPose = sva::PTransformd(Eigen::Vector3d{0.0, 0.0, 0.0}));

void sortFaceVertices(std::vector<Eigen::Vector3d> & vertices, Eigen::Vector3d faceNormal);

} // namespace mc_dynamic_polytopes
