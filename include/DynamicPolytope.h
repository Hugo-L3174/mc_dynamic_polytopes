#pragma once

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>

#include <mc_control/fsm/Controller.h>
#include <mc_rtc/clock.h>
#include <mc_rtc/gui.h>

#include "WrenchCones.h"

#include <politopix/PolyhedralAlgorithms_Rn.h>
#include <politopix/PrismaticPolyhedron_Rn.h>
#include <politopix/politopixAPI.h>

// #include <eigen-cdd/Polyhedron.h>

struct DynamicPolytope
{
  DynamicPolytope(const std::string & name, std::set<std::string> contactNames);
  ~DynamicPolytope();

  // compute the contact force cone into a 3D polytope
  void buildForceConeFromContact(int numberOfFrictionSides,
                                 std::pair<std::pair<double, double>, sva::PTransformd> & contactSurface,
                                 boost::shared_ptr<Polytope_Rn> & forceCone,
                                 double m_frictionCoef,
                                 double maxForce);

  // compute the force cone and the associated moment cone separately as two 3D polytopes
  void buildWrenchConeFromContact(int numberOfFrictionSides,
                                  std::pair<std::pair<double, double>, sva::PTransformd> & contactSurface,
                                  boost::shared_ptr<Polytope_Rn> & forceCone,
                                  boost::shared_ptr<Polytope_Rn> & momentPoly,
                                  double m_frictionCoef,
                                  double maxForce,
                                  Eigen::Vector3d CoM);

  // compute the force polytope of the contact from the wrench limits of the limb actuating it
  // args should be a polytope pointer/ref and a way to structure the limb Jacobian (only for a chosen limb)
  // think about how to compute bounds if jac non diagonal (redundancy): orso2018ral says QP? or LP?
  // also if not considering zero vel and acc how to use inertia matrix elements
  void buildActuationPolytopeFromContact(boost::shared_ptr<Polytope_Rn> & actuationPolytope);

  /* computes all cones from the surfaces with the given names (set by setCurrentContacts), reset the pointers of the
  map and updates the H-description of the poly using the double description algorithm.
  */
  void computeConesFromContactSet(const mc_rbdyn::Robot & robot);

  /* computes directly the V-rep of the CWC from individual contact friction cones and the moment limits transformed to
  the CoM (transformation from contact) then runs double description to update H-rep
  */
  void computeCWCFromContactSet(const mc_rbdyn::Robot & robot);

  /* Computes the 3d volume formed between the possible ZMP area(s) and the CoM of the robot
  TODO this is potentially several convex areas! (caron tro) see how to handle this
  */
  void computeZMPRegion(Eigen::Vector3d comPosition, const mc_rbdyn::Robot & robot);

  // Creates a 6d contact friction cone from the contact surface border points
  // The generators are computed then used to build the Polytope_Rn object which is added to the cones vector
  // void computeWrenchConesFromContactSet(const mc_rbdyn::Robot & robot);

  // Computes the minkowsky sum of the given friction cones and puts the result in the CWC_ polytope object
  void computeMinkowskySumPolitopix();

  void computeECMPRegion(Eigen::Vector3d comPosition, const mc_rbdyn::Robot & robot);

  void computeMomentsRegion(Eigen::Vector3d comPosition, const mc_rbdyn::Robot & robot);

  /* Computes the intersection between the eCMP region and the ZMP region to get the zero moment region,
  put in zeroMomentRegion_
  */
  void computeZeroMomentIntersection();

  // Computes the convex hull of the CWC_ polytope
  // Might be unnecessary, heavy algorithm to remove unnecessary faces
  void computeResultHull();

  // Updates the internal maps of triangles for gui display for the given contact names
  void updateTrianglesGUIPolitopix();

  std::vector<std::array<Eigen::Vector3d, 3>> getForceConesTriangles(const std::string & name)
  {
    return forceConesTrianglesMap_.at(name);
  };

  std::vector<std::array<Eigen::Vector3d, 3>> getContactMomentTriangles(const std::string & name)
  {
    return momentPolytopesTrianglesMap_.at(name);
  };

  std::vector<std::array<Eigen::Vector3d, 3>> getCWCForceTriangles()
  {
    return CWCForceTriangles_;
  };

  std::vector<std::array<Eigen::Vector3d, 3>> getCWCMomentTriangles()
  {
    return CWCMomentTriangles_;
  };

  std::vector<std::array<Eigen::Vector3d, 3>> getECMPTriangles()
  {
    return eCMPTriangles_;
  };

  std::vector<std::array<Eigen::Vector3d, 3>> getZMPTriangles()
  {
    return ZMPTriangles_;
  };

  std::vector<std::array<Eigen::Vector3d, 3>> getZeroMomentTriangles()
  {
    return zeroMomentTriangles_;
  };

  // From the current contact set, deduce what contacts need to be removed from computation compared to last iteration
  void setCurrentContacts(std::vector<std::string> & contactNames)
  {
    // take the previous set of contacts
    contactsToRemove_ = activeContacts_;
    activeContacts_.clear();
    for(const auto contactName : contactNames)
    {
      // add every contact given to the active contacts
      activeContacts_.emplace(contactName);
      // every active contact does not need to be removed
      contactsToRemove_.erase(contactName);
    }
  };

  void setWithMoments(bool withMoments)
  {
    withMoments_ = withMoments;
  };

  void computeECMP(const mc_rbdyn::Robot & robot);

  void load(const mc_rtc::Configuration & config);
  void addToGUI(mc_rtc::gui::StateBuilder & gui,
                double guiScale,
                std::vector<std::string> category = {"DynamicPolytopes"});
  void removeFromGUI(mc_rtc::gui::StateBuilder & gui);
  void addToLogger(mc_rtc::Logger & logger, const std::string & prefix = "DynamicPolytopes");
  void removeFromLogger(mc_rtc::Logger & logger);

protected:
  // Updates the faces vector used for polytope display (internal function)
  void update3DPolyTrianglesPolitopix(boost::shared_ptr<Polytope_Rn> & polytope,
                                      std::vector<std::array<Eigen::Vector3d, 3>> & resultTriangles,
                                      double guiScale);

  void update6DPolyTrianglesPolitopix(boost::shared_ptr<Polytope_Rn> & polytope,
                                      std::vector<std::array<Eigen::Vector3d, 3>> & resultMomentTriangles,
                                      std::vector<std::array<Eigen::Vector3d, 3>> & resultForceTriangles,
                                      double guiScale);

  void clearTriangles(std::vector<std::array<Eigen::Vector3d, 3>> & resultTriangles)
  {
    resultTriangles.clear();
  };

  std::string name_;
  std::set<std::string> possibleContacts_;
  std::set<std::string> activeContacts_;
  std::set<std::string> contactsToRemove_;

  bool withMoments_;

  mc_rtc::gui::PolyhedronConfig polyForceConfig_;
  mc_rtc::gui::PolyhedronConfig polyMomentConfig_;
  mc_rtc::gui::PolyhedronConfig polyZMPConfig_;
  mc_rtc::gui::PolyhedronConfig polyZeroMomentAreaConfig_;
  double guiScale_;

  sva::ForceVecd robotNetWrench_;
  Eigen::Vector3d eCMP_;

  // politopix
  std::map<std::string, boost::shared_ptr<Polytope_Rn>> frictionCones_;
  std::map<std::string, boost::shared_ptr<Polytope_Rn>> frictionConesMoments_;

  boost::shared_ptr<Polytope_Rn> CWCForces_;
  boost::shared_ptr<Polytope_Rn> CWCMoments_;
  // boost::shared_ptr<Polytope_Rn> eCMPRegion_;

  boost::shared_ptr<Polytope_Rn> zmpRegion_;
  boost::shared_ptr<Polytope_Rn> zeroMomentRegion_;

  // cdd
  // std::vector<std::shared_ptr<Eigen::Polyhedron>> cddFrictionCones_;

  // map of polytope triangles for display
  std::map<std::string, std::vector<std::array<Eigen::Vector3d, 3>>> forceConesTrianglesMap_;
  std::map<std::string, std::vector<std::array<Eigen::Vector3d, 3>>> momentPolytopesTrianglesMap_;
  std::vector<std::array<Eigen::Vector3d, 3>> CWCForceTriangles_;
  std::vector<std::array<Eigen::Vector3d, 3>> CWCMomentTriangles_;
  std::vector<std::array<Eigen::Vector3d, 3>> eCMPTriangles_;
  std::vector<std::array<Eigen::Vector3d, 3>> ZMPTriangles_;
  std::vector<std::array<Eigen::Vector3d, 3>> zeroMomentTriangles_;
};
