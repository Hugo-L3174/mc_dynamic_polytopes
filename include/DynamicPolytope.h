#pragma once

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>

#include <mc_control/fsm/Controller.h>
#include <mc_rtc/clock.h>
#include <mc_rtc/gui.h>

#include "GUIComputations.h"
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

  // Set current contact set to be used for next computation
  void setControllerContacts(const std::vector<std::string> & contactNames)
  {
    auto waitForLock = mc_rtc::clock::now();
    std::lock_guard<std::mutex> lock(contactSetMutex_);
    mc_rtc::duration_ms lockTime = mc_rtc::clock::now() - waitForLock;
    // mc_rtc::log::critical("Controller waited {}ms to get contact lock", lockTime.count());
    controllerContactsNames_ = contactNames;
    // if computation has not started yet, launch it
    if(!computing_)
    {
      computing_ = true;
      // signal main computation thread
      cv_.notify_one();
    }
  }

  // Set current contact set to be used for next computation
  void setControllerContacts(const std::vector<mc_rbdyn::Contact> & contacts)
  {
    auto waitForLock = mc_rtc::clock::now();
    std::lock_guard<std::mutex> lock(contactSetMutex_);
    mc_rtc::duration_ms lockTime = mc_rtc::clock::now() - waitForLock;
    // mc_rtc::log::critical("Controller waited {}ms to get contact lock", lockTime.count());
    controllerContacts_ = contacts;
    // if computation has not started yet, launch it
    if(!computing_)
    {
      computing_ = true;
      // signal main computation thread
      cv_.notify_one();
    }
  }

  // Set the use of moments for computations
  void setWithMoments(bool withMoments)
  {
    withMoments_ = withMoments;
  };

  Eigen::Vector3d computeECMP(const mc_rbdyn::Robot & robot);

  // ------------------------------------------------------> mc_rtc interface functions

  void load(const mc_rtc::Configuration & config);
  void addToGUI(mc_rtc::gui::StateBuilder & gui,
                double guiScale,
                std::vector<std::string> category = {"DynamicPolytopes"});
  void removeFromGUI(mc_rtc::gui::StateBuilder & gui);
  void addToLogger(mc_rtc::Logger & logger, const std::string & prefix = "DynamicPolytopes_");
  void removeFromLogger(mc_rtc::Logger & logger);

  // ------------------------------------------------------> public computation functions

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

  // ------------------------------------------------------> Plane constraints getters

  // Returns the internal normals matrix and offsets vector of the eCMP region for QP constraint or check
  void getECMPPlanes(Eigen::MatrixX3d & Normals, Eigen::VectorXd & Offsets)
  {
    std::lock_guard<std::mutex> lock(eCMPPlanesMutex_);
    Normals.resize(eCMPNormals_.rows(), 3);
    Offsets.resize(eCMPOffsets_.size());
    Normals = eCMPNormals_;
    Offsets = eCMPOffsets_;
  };

  // Returns the internal normals matrix and offsets vector of the zero moment region (subset of the DCM region) for QP
  // constraint or check
  void getZeroMomentPlanes(Eigen::MatrixX3d & Normals, Eigen::VectorXd & Offsets)
  {
    std::lock_guard<std::mutex> lock(zeroMomentPlanesMutex_);
    Normals.resize(zeroMomentNormals_.rows(), 3);
    Offsets.resize(zeroMomentOffsets_.size());
    Normals = zeroMomentNormals_;
    Offsets = zeroMomentOffsets_;
  };

  // Projects the given point on the VRP region. Returns the given point if it is already inside
  Eigen::Vector3d projectPointInVRPRegion(Eigen::Vector3d testedPoint);

  Eigen::Vector3d projectPointInZeroMomentRegion(Eigen::Vector3d testedPoint);

protected:
  // main computation function that calls all region calculations in sequence
  // intended to be called in a thread, and will thread region calculations correctly
  void computeRegions();

  // Updates the internal maps of triangles for gui display for the given contact names
  void updateTrianglesGUIPolitopix();

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
  void computeZMPRegion(const mc_rbdyn::Robot & robot,const std::vector<mc_rbdyn::Contact> & contacts);

  void computeMomentLessForceCone(const mc_rbdyn::Robot & robot,const std::vector<mc_rbdyn::Contact> & contacts);

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

  // Translates the eCMP region and the zero-moment region with a vertical Delta z offset to get actual DCM region
  void VRPtranslation(double deltaZ);

  // Computes the convex hull of the CWC_ polytope
  // Might be unnecessary, heavy algorithm to remove unnecessary faces
  void computeResultHull();

  // Puts the H representation of the given polytope into a matrix (normals) and a vector (offsets) for easy
  // testing/constraining. /!\ politopix convention has normals towards the inside, so we negate them again to return
  // them in usual convention (normals towards exterior)
  void updatePlanesMatrixConstraint(const boost::shared_ptr<Polytope_Rn> & polytope,
                                    Eigen::MatrixX3d & Normals,
                                    Eigen::VectorXd & Offsets);


  Eigen::Vector3d projectPointInPolytope(Eigen::Vector3d testedPoint, boost::shared_ptr<Polytope_Rn> & polytope);

  // sanity check
  bool checkGravityCenterInPolytope(boost::shared_ptr<Polytope_Rn> & polytope);

  // ------------------------------------------------------> Timing functions

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

  inline mc_rtc::duration_ms dt_zeroMomentInter() const noexcept
  {
    return dt_zeroMoment_intersection_;
  }

  // ------------------------------------------------------> GUI getters

  std::vector<std::array<Eigen::Vector3d, 3>> getForceConesTriangles(const std::string & name)
  {
    std::lock_guard<std::mutex> lock(getContactMutex(forceConeTrianglesMutexes_, name));
    return forceConesTrianglesMap_.at(name);
  };

  std::vector<std::array<Eigen::Vector3d, 3>> getContactMomentTriangles(const std::string & name)
  {
    std::lock_guard<std::mutex> lock(getContactMutex(momentTrianglesMutexes_, name));
    return momentPolytopesTrianglesMap_.at(name);
  };

  std::vector<std::array<Eigen::Vector3d, 3>> getCWCForceTriangles()
  {
    std::lock_guard<std::mutex> lock(CWCForceTrianglesMutex_);
    return CWCForceTriangles_;
  };

  std::vector<std::array<Eigen::Vector3d, 3>> getCWCMomentTriangles()
  {
    std::lock_guard<std::mutex> lock(CWCMomentTrianglesMutex_);
    return CWCMomentTriangles_;
  };

  std::vector<std::array<Eigen::Vector3d, 3>> getCWCMomentLess()
  {
    return CWCMomentLessUnCstrTriangles_;
  };

  
  // std::vector<std::array<Eigen::Vector3d, 3>> getECMPTriangles()
  // {
  //   return eCMPTriangles_;
  // };

  std::vector<std::array<Eigen::Vector3d, 3>> getZMPTriangles()
  {
    std::lock_guard<std::mutex> lock(ZMPTrianglesMutex_);
    return ZMPTriangles_;
  };

  std::vector<std::array<Eigen::Vector3d, 3>> getZeroMomentTriangles()
  {
    std::lock_guard<std::mutex> lock(zeroMomentTrianglesMutex_);
    return zeroMomentTriangles_;
  };

protected:

  // From the current contact set, deduce what contacts need to be removed from computation compared to last iteration
  void setCurrentContactsName()
  {
    auto waitForLock = mc_rtc::clock::now();
    std::lock_guard<std::mutex> lock(contactSetMutex_);
    mc_rtc::duration_ms lockTime = mc_rtc::clock::now() - waitForLock;
    // mc_rtc::log::critical("Region compute thread waited {}ms to get contact lock", lockTime.count());
    // take the previous set of contacts
    contactsToRemove_ = activeContacts_;
    activeContacts_.clear();
    for(const auto contactName : controllerContactsNames_)
    {
      // add every contact given to the active contacts
      activeContacts_.emplace(contactName);
      // every active contact does not need to be removed
      contactsToRemove_.erase(contactName);
    }
  };

  void setCurrentContacts()
  {
    auto waitForLock = mc_rtc::clock::now();
    std::lock_guard<std::mutex> lock(contactSetMutex_);
    mc_rtc::duration_ms lockTime = mc_rtc::clock::now() - waitForLock;
    // mc_rtc::log::critical("Region compute thread waited {}ms to get contact lock", lockTime.count());
    // take the previous set of contacts
    contacts_ = controllerContacts_; 

  };

  std::string name_;
  const mc_rbdyn::Robot & robot_;
  std::vector<mc_rbdyn::Contact> contacts_;
  std::vector<mc_rbdyn::Contact> controllerContacts_;
  std::set<std::string> possibleContacts_;
  std::set<std::string> activeContacts_;
  std::set<std::string> contactsToRemove_;

  // intermediate names vector to be set by controller and used once per compute loop
  std::vector<std::string> controllerContactsNames_;

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
  boost::shared_ptr<Polytope_Rn> CWCMomentLessUnCstr_;

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
  std::thread zmpThread_;
  std::mutex eCMPPlanesMutex_;
  std::mutex zeroMomentPlanesMutex_;
  // this is a workaround for the fact that mutexes are not movable
  // using this function to get the desired mutex from the name they will be created when needed
  std::mutex & getContactMutex(std::map<std::string, std::mutex> & mutexMap, const std::string & contactName)
  {
    return mutexMap[contactName]; // constructs it inside the map if doesn't exist
  }
  // gui mutexes
  std::map<std::string, std::mutex> forceConeTrianglesMutexes_;
  std::map<std::string, std::mutex> momentTrianglesMutexes_;
  std::mutex CWCForceTrianglesMutex_;
  std::mutex CWCMomentTrianglesMutex_;
  std::mutex ZMPTrianglesMutex_;
  std::mutex zeroMomentTrianglesMutex_;
  std::mutex CWCMomentLessUnCstrTrianglesMutex_;

  // Internal matrices of planes and offsets of the regions for constraints
  Eigen::MatrixX3d eCMPNormals_;
  Eigen::VectorXd eCMPOffsets_;

  Eigen::MatrixX3d zeroMomentNormals_;
  Eigen::VectorXd zeroMomentOffsets_;

  // timers to measure computation times
  mc_rtc::duration_ms dt_loop_total_;
  mc_rtc::duration_ms dt_compute_contactSet_;
  mc_rtc::duration_ms dt_compute_minkSum_;
  mc_rtc::duration_ms dt_update_planes_;
  mc_rtc::duration_ms dt_compute_guiTriangles_;
  mc_rtc::duration_ms dt_zeroMoment_intersection_;

  // map of polytope triangles for display
  std::map<std::string, std::vector<std::array<Eigen::Vector3d, 3>>> forceConesTrianglesMap_;
  std::map<std::string, std::vector<std::array<Eigen::Vector3d, 3>>> momentPolytopesTrianglesMap_;
  std::vector<std::array<Eigen::Vector3d, 3>> CWCForceTriangles_;
  std::vector<std::array<Eigen::Vector3d, 3>> CWCMomentTriangles_;
  std::vector<std::array<Eigen::Vector3d, 3>> CWCMomentLessUnCstrTriangles_;
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
