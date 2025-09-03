#pragma once
#include <SpaceVecAlg/SpaceVecAlg>
#include <politopix/politopixAPI.h>

// TODO: come up with a better filename and/or organization

namespace mc_dynamic_polytopes
{

// Computes the convex hull of the set of points given (Qhull format) and builds the halfspaces in the given polytope
void computeQhullHrep(std::vector<double> & points, boost::shared_ptr<Polytope_Rn> & polytope, int dim);

// compute the contact friction cone into an unbounded polyhedral cone (only planes in a polytope object)
void buildFrictionConeFromContactWithHrep(int numberOfFrictionSides,
                                          const sva::PTransformd X_r1_r2,
                                          boost::shared_ptr<Polytope_Rn> & frictionCone,
                                          double m_frictionCoef,
                                          unsigned int dim);
// Puts the H representation of the given polytope into a matrix (normals) and a vector (offsets) for easy
// testing/constraining. /!\ politopix convention has normals towards the inside, so we negate them again to return
// them in usual convention (normals towards exterior)
// TODO template this for polyhedral cones
void updatePlanesMatrixConstraint(const boost::shared_ptr<Polytope_Rn> & polytope,
                                  Eigen::MatrixXd & Normals,
                                  Eigen::VectorXd & Offsets);

void checkAllHSInternal(const std::string & polyName, boost::shared_ptr<Polytope_Rn> & polytope);

} // namespace mc_dynamic_polytopes
