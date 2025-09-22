#pragma once

#include <mc_rtc/gui.h>
#include <politopix/politopixAPI.h>
#include <vector>

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

/**
 * @brief Sorts the vertices of a face in counter-clockwise order around their centroid, projected onto the 2D face
 * plane.
 *
 * For faces with more than three vertices, this function sorts the input vertices in-place so that they are ordered
 * counter-clockwise around the centroid of the face, as seen along the provided face normal. The centroid is computed
 * and returned. For triangular faces (three vertices), the function ensures the vertices are ordered counter-clockwise
 * with respect to the face normal. If there are fewer than three vertices, the function returns std::nullopt and does
 * not modify the input.
 *
 * @param vertices Vector of 3D points representing the face vertices (modified in-place).
 * @param faceNormal The normal vector of the face (should be normalized).
 * @return The centroid of the face if sorting was performed (for n > 3), std::nullopt otherwise.
 */
std::optional<Eigen::Vector3d> sortFaceVertices(std::vector<Eigen::Vector3d> & vertices, Eigen::Vector3d faceNormal);

} // namespace mc_dynamic_polytopes
