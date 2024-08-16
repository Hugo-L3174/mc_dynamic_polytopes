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
  DynamicPolytope(const std::string & name, std::set<std::string> contactNames, const mc_rbdyn::Robot & robot);
  ~DynamicPolytope();

  // stopping function for destructor
  inline void stopThread()
  {
    // XXX stop inner threads as well?
    computing_ = false;
    mainComputeThread_.join();
  }

  // main computation function that calls all region calculations in sequence
  // intended to be called in a thread, and will thread region calculations correctly
  void computeRegions();

  // compute the contact force cone into a 3D polytope
  void buildForceConeFromContact(int numberOfFrictionSides,
                                 sva::PTransformd contactSurface,
                                 boost::shared_ptr<Polytope_Rn> & forceCone,
                                 std::mutex & forceConeMutex,
                                 double m_frictionCoef,
                                 double maxForce);

  // compute the force cone and the associated moment cone separately as two 3D polytopes
  void buildWrenchConeFromContact(int numberOfFrictionSides,
                                  std::pair<std::pair<double, double>, sva::PTransformd> contactSurface,
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
  void computeZMPRegion(Eigen::Vector3d comPosition);

  // Creates a 6d contact friction cone from the contact surface border points
  // The generators are computed then used to build the Polytope_Rn object which is added to the cones vector
  // void computeWrenchConesFromContactSet(const mc_rbdyn::Robot & robot);

  // Computes the minkowsky sum of the given friction cones and puts the result in the CWC_ polytope object
  void computeMinkowskySumPolitopix();

  // Computes the eCMP region from the CWC
  // /!\ this modifies the CWC to invert and scale it, the result is still in the CWC polytope !
  void computeECMPRegion(Eigen::Vector3d comPosition, const mc_rbdyn::Robot & robot);

  // Computes the 3d volume of the moments from the contact surfaces borders resulting from the CWC generating rays at
  // these limits
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

  // Updates the given normals matrix and offsets vector of the eCMP region for QP constraint or check
  void updateECMPPlanes(Eigen::MatrixX3d & Normals, Eigen::VectorXd & Offsets)
  {
    updatePlanesMatrixConstraint(CWCForces_, Normals, Offsets);
  };

  // Updates the given normals matrix and offsets vector of the zero moment region (subset of the DCM region) for QP
  // constraint or check
  void updateZeroMomentPlanes(Eigen::MatrixX3d & Normals, Eigen::VectorXd & Offsets)
  {
    updatePlanesMatrixConstraint(zeroMomentRegion_, Normals, Offsets);
  };

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

  // Set current contact set to be used for next computation
  void setControllerContacts(const std::vector<std::string> & contactNames)
  {
    auto waitForLock = mc_rtc::clock::now();
    std::lock_guard<std::mutex> lock(contactSetMutex_);
    mc_rtc::duration_ms lockTime = mc_rtc::clock::now() - waitForLock;
    mc_rtc::log::critical("Controller waited {}ms to get contact lock", lockTime.count());
    controllerContacts_ = contactNames;
    // if computation has not started yet, launch it
    if(!computing_)
    {
      computing_ = true;
      // signal main computation thread
      cv_.notify_one();
    }
  }

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
  void addToLogger(mc_rtc::Logger & logger, const std::string & prefix = "DynamicPolytopes_");
  void removeFromLogger(mc_rtc::Logger & logger);

  inline mc_rtc::duration_ms dt_loop_total() const noexcept
  {
    return dt_loop_total_;
  }

  inline mc_rtc::duration_ms dt_contactSet() const noexcept
  {
    return dt_compute_contactSet_;
  }

  inline mc_rtc::duration_ms dt_minkSum() const noexcept
  {
    return dt_compute_minkSum_;
  }

  inline mc_rtc::duration_ms dt_updatePlanes() const noexcept
  {
    return dt_update_planes_;
  }

  inline mc_rtc::duration_ms dt_guiTriangles() const noexcept
  {
    return dt_compute_guiTriangles_;
  }

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

  // Puts the H representation of the given polytope into a matrix (normals) and a vector (offsets) for easy
  // testing/constraining
  void updatePlanesMatrixConstraint(const boost::shared_ptr<Polytope_Rn> & polytope,
                                    Eigen::MatrixX3d & Normals,
                                    Eigen::VectorXd & Offsets);

  // From the current contact set, deduce what contacts need to be removed from computation compared to last iteration
  void setCurrentContacts()
  {
    auto waitForLock = mc_rtc::clock::now();
    std::lock_guard<std::mutex> lock(contactSetMutex_);
    mc_rtc::duration_ms lockTime = mc_rtc::clock::now() - waitForLock;
    mc_rtc::log::critical("Region compute thread waited {}ms to get contact lock", lockTime.count());
    // take the previous set of contacts
    contactsToRemove_ = activeContacts_;
    activeContacts_.clear();
    for(const auto contactName : controllerContacts_)
    {
      // add every contact given to the active contacts
      activeContacts_.emplace(contactName);
      // every active contact does not need to be removed
      contactsToRemove_.erase(contactName);
    }
  };

  std::string name_;
  const mc_rbdyn::Robot & robot_;
  std::set<std::string> possibleContacts_;
  std::set<std::string> activeContacts_;
  std::set<std::string> contactsToRemove_;

  // intermediate names vector to be set by controller and used once per compute loop
  std::vector<std::string> controllerContacts_;

  bool withMoments_;

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

  // threading things
  std::thread mainComputeThread_;
  // atomic bool as condition to keep computing in loop
  std::atomic<bool> computing_{false};
  // condition variable to signal start of loop for main thread
  std::condition_variable cv_;
  std::mutex contactSetMutex_;
  std::map<std::string, std::thread> frictionConesThreads_;
  std::map<std::string, std::mutex> frictionConesMutexes_;
  // this is a workaround for the fact that mutexes are not movable
  // using this function to get the desired mutex from the name they will be created when needed
  std::mutex & getContactMutex(const std::string & contactName)
  {
    return frictionConesMutexes_[contactName]; // constructs it inside the map if doesn't exist
  }

  // timers to measure computation times
  mc_rtc::duration_ms dt_loop_total_;
  mc_rtc::duration_ms dt_compute_contactSet_;
  mc_rtc::duration_ms dt_compute_minkSum_;
  mc_rtc::duration_ms dt_update_planes_;
  mc_rtc::duration_ms dt_compute_guiTriangles_;

  // map of polytope triangles for display
  std::map<std::string, std::vector<std::array<Eigen::Vector3d, 3>>> forceConesTrianglesMap_;
  std::map<std::string, std::vector<std::array<Eigen::Vector3d, 3>>> momentPolytopesTrianglesMap_;
  std::vector<std::array<Eigen::Vector3d, 3>> CWCForceTriangles_;
  std::vector<std::array<Eigen::Vector3d, 3>> CWCMomentTriangles_;
  std::vector<std::array<Eigen::Vector3d, 3>> eCMPTriangles_;
  std::vector<std::array<Eigen::Vector3d, 3>> ZMPTriangles_;
  std::vector<std::array<Eigen::Vector3d, 3>> zeroMomentTriangles_;

  // gui configs
  mc_rtc::gui::PolyhedronConfig polyForceConfig_;
  mc_rtc::gui::PolyhedronConfig polyMomentConfig_;
  mc_rtc::gui::PolyhedronConfig polyZMPConfig_;
  mc_rtc::gui::PolyhedronConfig polyZeroMomentAreaConfig_;
  double guiScale_;
};
