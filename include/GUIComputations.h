#pragma once

#include <mc_rtc/gui.h>
#include <politopix/politopixAPI.h>


// Updates the faces vector used for polytope display
void update3DPolyTrianglesPolitopix(boost::shared_ptr<Polytope_Rn> & polytope,
                                    std::vector<std::array<Eigen::Vector3d, 3>> & resultTriangles,
                                    double guiScale);

// Updates the faces vector of the 3D force and 3D moment in the case of a 6D polytope
void update6DPolyTrianglesPolitopix(boost::shared_ptr<Polytope_Rn> & polytope,
                                    std::vector<std::array<Eigen::Vector3d, 3>> & resultMomentTriangles,
                                    std::vector<std::array<Eigen::Vector3d, 3>> & resultForceTriangles,
                                    double guiScale);
