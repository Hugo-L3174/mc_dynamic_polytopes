#pragma once

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>

#include <RBDyn/Coriolis.h>

#include <mc_control/fsm/Controller.h>
#include <mc_rtc/clock.h>
#include <mc_rtc/gui.h>

#include "GUIComputations.h"
#include "WrenchCones.h"

#include <politopix/PolyhedralAlgorithms_Rn.h>
#include <politopix/PrismaticPolyhedron_Rn.h>
#include <politopix/politopixAPI.h>

#include <libqhull_r/qhull_ra.h>

using HRepX3d = std::pair<Eigen::MatrixX3d, Eigen::VectorXd>;

struct ContactTimers
{
  mc_rtc::duration_ms dt_frictionCone;
  mc_rtc::duration_ms dt_forcePolytope;
  mc_rtc::duration_ms dt_intersection;
  mc_rtc::duration_ms dt_contactTotal;
};

struct DynamicPolytope
{
  DynamicPolytope(const std::string & name,
                  const mc_rbdyn::Robot & robot,
                  const mc_rtc::Configuration & dynamicPolyConfig);
  ~DynamicPolytope();

  // stopping function for destructor
  inline void stopThread()
  {
    computing_ = false;
    mainComputeThread_.join();
    if(zmpThread_.joinable())
    {
      zmpThread_.join();
    }
    if(minkSumThread_.joinable())
    {
      minkSumThread_.join();
    }
    for(auto & thread : feasiblePolytopesThreads_)
    {
      thread.second.join();
    }
  }

  // Set current contact set to be used for next computation: contact name and contact reference pose, ie the frame of
  // the target surface of the contact (only the orientation is used for the friction cone)
  void setControllerContacts(const std::vector<std::pair<std::string, sva::PTransformd>> & contacts)
  {
    auto waitForLock = mc_rtc::clock::now();
    std::lock_guard<std::mutex> lock(contactSetMutex_);
    mc_rtc::duration_ms lockTime = mc_rtc::clock::now() - waitForLock;
    // mc_rtc::log::critical("Controller waited {}ms to get contact lock", lockTime.count());
    controllerContacts_.clear();
    refContactTransforms_.clear();
    for(const auto & contact : contacts)
    {
      controllerContacts_.emplace_back(contact.first);
      refContactTransforms_.emplace(contact.first, contact.second);
    }

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

  // ------------------------------------------------------> mc_rtc interface functions

  void addToGUI(mc_rtc::gui::StateBuilder & gui,
                double guiScale = 0.001,
                std::vector<std::string> category = {"DynamicPolytopes"});
  void removeFromGUI(mc_rtc::gui::StateBuilder & gui);
  void addToLogger(mc_rtc::Logger & logger, const std::string & prefix = "DynamicPolytopes_");
  void removeFromLogger(mc_rtc::Logger & logger);

  // ------------------------------------------------------> public computation functions

  // compute the contact friction cone into an unbounded polyhedral cone (only planes in a polytope object)
  void buildFrictionConeFromContactWithHrep(int numberOfFrictionSides,
                                            const sva::PTransformd X_r1_r2,
                                            boost::shared_ptr<Polytope_Rn> & frictionCone,
                                            std::mutex & frictionConeMutex,
                                            double m_frictionCoef);

  // compute the contact friction cone into a bounded polytope from the Vrep, with planes formed from the bounding
  // volume computation
  void buildFrictionConeFromContactWithVrep(int numberOfFrictionSides,
                                            const sva::PTransformd contactSurface,
                                            boost::shared_ptr<Polytope_Rn> & frictionCone,
                                            std::mutex & frictionConeMutex,
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
  // now also taking a scale factor (0-1) for force
  void buildActuationPolytopeFromContact(const std::string contactName,
                                         const mc_rbdyn::Robot & robot,
                                         boost::shared_ptr<Polytope_Rn> & actuationPolytope,
                                         std::mutex & forcePolyMutex,
                                         double forceScalingFactor);

  // compute friction cone and force polytope, then intersect both into the friction cone object
  void buildFeasiblePolytopeFromContact(const std::string contactName,
                                        const mc_rbdyn::Robot & robot,
                                        const sva::PTransformd refContactPose,
                                        int numberOfFrictionSides,
                                        double forceScalingFactor,
                                        double m_frictionCoef,
                                        boost::shared_ptr<Polytope_Rn> & frictionCone,
                                        std::mutex & frictionConeMutex,
                                        boost::shared_ptr<Polytope_Rn> & actuationPolytope,
                                        std::mutex & forcePolyMutex,
                                        ContactTimers & timers);

  // ------------------------------------------------------> Plane constraints getters

  // Returns the internal normals matrix and offsets vector of the eCMP region for QP constraint or check
  const HRepX3d & getVRPPlanes()
  {
    std::lock_guard<std::mutex> lock(VRPPlanesMutex_);
    return DCMVRPPlanes_;
  };

  // Returns the internal normals matrix and offsets vector of the zero moment region (subset of the DCM region) for QP
  // constraint or check
  const HRepX3d & getZeroMomentPlanes()
  {
    std::lock_guard<std::mutex> lock(zeroMomentPlanesMutex_);
    return zeroMomentPlanes_;
  };

  // Returns the desired cone H representation
  // Prefer using the force polytope planes if they are computed
  const HRepX3d & getFrictionConePlanes(const std::string & contactName)
  {
    std::lock_guard<std::mutex> lock(getContactMutex(frictionConesPlanesMutexes_, contactName));
    return frictionConesPlanes_.at(contactName);
  };

  // Returns the desired contact force polytope H representation
  const HRepX3d & getForcePolyPlanes(const std::string & contactName)
  {
    std::lock_guard<std::mutex> lock(getContactMutex(forcePolyPlanesMutexes_, contactName));
    return forcePolyPlanes_.at(contactName);
  };

  // Projects the given point on the VRP region. Returns the given point if it is already inside
  Eigen::Vector3d projectPointInVRPRegion(Eigen::Vector3d testedPoint);

  Eigen::Vector3d projectPointInZeroMomentRegion(Eigen::Vector3d testedPoint);

protected:
  // main computation function that calls all region calculations in sequence
  // intended to be called in a thread, and will thread region calculations correctly
  void computeRegions();

  // Updates the internal maps of triangles of the contacts for gui display
  void updateTrianglesContactsGUIPolitopix();

  // Updates the internal maps of triangles of the total regions for gui display
  void updateTrianglesRegionsGUIPolitopix();

  /* computes all cones from the surfaces with the given names (set by setCurrentContacts), reset the pointers of the
  map and updates the H-description of the cone using the double description algorithm.
  */
  void computeFrictionConesFromContactSet(const mc_rbdyn::Robot & robot);

  /* computes all force polytopes from the surfaces with the given names (set by setCurrentContacts), reset the pointers
  of the map and updates the V-description of the poly using the double description algorithm.
  */
  void computeForcePolyFromContactSet(const mc_rbdyn::Robot & robot);

  /* Computes for each contact in parallel: the friction cone and force polytope,
  then their intersection back into the friction cone polytope
  */
  void computeFeasibleForcesFromContactSet(const mc_rbdyn::Robot & robot);

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

  // Translates the eCMP region and the zero-moment region with a vertical Delta z offset to get actual DCM region
  void VRPtranslation(double deltaZ);

  // Computes the mink sum of the contact regions, scale to get eCMP region, intersect with ZMP region, and translate
  // everything to get VRP region
  void computeVRPRegionWithMinkSum();

  // Computes the convex hull of the CWC_ polytope
  // Might be unnecessary, heavy algorithm to remove unnecessary faces
  void computeResultHull();

  // Puts the H representation of the given polytope into a matrix (normals) and a vector (offsets) for easy
  // testing/constraining. /!\ politopix convention has normals towards the inside, so we negate them again to return
  // them in usual convention (normals towards exterior)
  // TODO template this for polyhedral cones
  void updatePlanesMatrixConstraint(const boost::shared_ptr<Polytope_Rn> & polytope,
                                    Eigen::MatrixX3d & Normals,
                                    Eigen::VectorXd & Offsets);

  // From the current contact set, deduce what contacts need to be removed from computation compared to last iteration
  void setCurrentContacts()
  {
    auto waitForLock = mc_rtc::clock::now();
    std::lock_guard<std::mutex> lock(contactSetMutex_);
    mc_rtc::duration_ms lockTime = mc_rtc::clock::now() - waitForLock;
    // mc_rtc::log::critical("Region compute thread waited {}ms to get contact lock", lockTime.count());
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

  // Computes the convex hull of the set of points given (Qhull format) and builds the halfspaces in the given polytope
  void computeQhullHrep(std::vector<double> & points, boost::shared_ptr<Polytope_Rn> & polytope);

  Eigen::Vector3d projectPointInPolytope(Eigen::Vector3d testedPoint, boost::shared_ptr<Polytope_Rn> & polytope);

  // sanity check
  bool checkGravityCenterInPolytope(boost::shared_ptr<Polytope_Rn> & polytope);

  void checkAllHSInternal(const std::string & polyName, boost::shared_ptr<Polytope_Rn> & polytope);

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

  inline mc_rtc::duration_ms dt_guiTrianglesContacts() const noexcept
  {
    return dt_compute_guiTrianglesContacts_;
  }

  inline mc_rtc::duration_ms dt_guiTrianglesRegions() const noexcept
  {
    return dt_compute_guiTrianglesRegions_;
  }

  inline mc_rtc::duration_ms dt_zeroMomentInter() const noexcept
  {
    return dt_zeroMoment_intersection_;
  }

  inline mc_rtc::duration_ms getContact_dt_frictionCone(const std::string & name) const noexcept
  {
    return contactsTimers_.at(name).dt_frictionCone;
  }

  inline mc_rtc::duration_ms getContact_dt_forcePolytope(const std::string & name) const noexcept
  {
    return contactsTimers_.at(name).dt_forcePolytope;
  }

  inline mc_rtc::duration_ms getContact_dt_intersection(const std::string & name) const noexcept
  {
    return contactsTimers_.at(name).dt_intersection;
  }

  inline mc_rtc::duration_ms getContact_dt_Total(const std::string & name) const noexcept
  {
    return contactsTimers_.at(name).dt_contactTotal;
  }

  // ------------------------------------------------------> GUI getters

  std::vector<std::array<Eigen::Vector3d, 3>> getFrictionConesTriangles(const std::string & name)
  {
    std::lock_guard<std::mutex> lock(getContactMutex(frictionConeTrianglesMutexes_, name));
    return frictionConesTrianglesMap_.at(name);
  };

  std::vector<std::array<Eigen::Vector3d, 3>> getForcePolyTriangles(const std::string & name)
  {
    std::lock_guard<std::mutex> lock(getContactMutex(forcePolyTrianglesMutexes_, name));
    return forcePolyTrianglesMap_.at(name);
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

  double & getForceScalingFactor(const std::string & name)
  {
    return forceScalingFactors_.at(name);
  }

  double & getFrictionCoeff(const std::string & name)
  {
    return frictionCoefficients_.at(name);
  }

  int & getNbOfFrictionSides(const std::string & name)
  {
    return numbersOfFrictionSides_.at(name);
  }

  // ------------------------------------------------------> Internal variables

  mc_rtc::Configuration config_;
  std::string name_;
  const mc_rbdyn::Robot & robot_;
  std::set<std::string> possibleContacts_;
  std::set<std::string> activeContacts_;
  std::set<std::string> contactsToRemove_;

  // map of the contact transforms X_r1_r2, with r1 the controlled robot frame and r2 the target frame of the contact
  // This is used to compute the orientation of the friction cone as it should be in the r1 frame (same as the force
  // polytope) for distribution, and not in world frame
  std::map<std::string, sva::PTransformd> refContactTransforms_;

  // map of the force scaling factors (alphas to be used to transfer between contacts)
  std::map<std::string, double> forceScalingFactors_;

  // friction coeffs for the cones, stored individually from contact name
  std::map<std::string, double> frictionCoefficients_;
  // number of sides for cones linearization, stored individually from contact name
  std::map<std::string, int> numbersOfFrictionSides_;

  // intermediate names vector to be set by controller and used once per compute loop
  std::vector<std::string> controllerContacts_;

  bool withMoments_;
  bool computeRegions_;
  bool HrepMode_;

  // politopix

  // We keep the friction cones as polytope objects, but they are actually polyhedral cones and will not be bounded
  std::map<std::string, boost::shared_ptr<Polytope_Rn>> frictionCones_;
  // The bounded actuation polytopes
  std::map<std::string, boost::shared_ptr<Polytope_Rn>> forcePolytopes_;
  std::map<std::string, boost::shared_ptr<Polytope_Rn>> frictionConesMoments_;

  boost::shared_ptr<Polytope_Rn> CWCForces_;
  boost::shared_ptr<Polytope_Rn> CWCMoments_;
  // boost::shared_ptr<Polytope_Rn> eCMPRegion_;

  boost::shared_ptr<Polytope_Rn> zmpRegion_;
  boost::shared_ptr<Polytope_Rn> zeroMomentRegion_;

  // threading things
  std::thread mainComputeThread_;
  // atomic bool as condition to keep computing in loop
  std::atomic<bool> computing_{false};
  // condition variable to signal start of loop for main thread
  std::condition_variable cv_;
  std::mutex contactSetMutex_;
  std::map<std::string, std::thread> feasiblePolytopesThreads_;
  std::mutex feasiblePolytopesThreadsMutex_;
  std::map<std::string, std::mutex> frictionConesMutexes_;
  std::map<std::string, std::mutex> forcePolyMutexes_;
  std::thread zmpThread_;
  std::thread minkSumThread_;
  std::atomic<bool> VRPRegionComputed_{true};
  std::mutex VRPPlanesMutex_;
  std::mutex zeroMomentPlanesMutex_;
  std::map<std::string, std::mutex> frictionConesPlanesMutexes_;
  std::map<std::string, std::mutex> forcePolyPlanesMutexes_;

  std::mutex CWCMutex_;
  std::mutex ZMPMutex_;
  std::mutex zeroMomentMutex_;

  // this is a workaround for the fact that mutexes are not movable
  // using this function to get the desired mutex from the name they will be created when needed
  std::mutex & getContactMutex(std::map<std::string, std::mutex> & mutexMap, const std::string & contactName)
  {
    return mutexMap[contactName]; // constructs it inside the map if doesn't exist
  }
  // gui mutexes
  std::map<std::string, std::mutex> frictionConeTrianglesMutexes_;
  std::map<std::string, std::mutex> forcePolyTrianglesMutexes_;
  std::map<std::string, std::mutex> momentTrianglesMutexes_;
  std::mutex CWCForceTrianglesMutex_;
  std::mutex CWCMomentTrianglesMutex_;
  std::mutex ZMPTrianglesMutex_;
  std::mutex zeroMomentTrianglesMutex_;

  // Internal matrices of planes and offsets of the regions for constraints
  HRepX3d DCMVRPPlanes_;

  HRepX3d zeroMomentPlanes_;

  std::map<std::string, HRepX3d> frictionConesPlanes_;
  std::map<std::string, HRepX3d> forcePolyPlanes_;

  // timers to measure computation times
  mc_rtc::duration_ms dt_loop_total_;
  mc_rtc::duration_ms dt_compute_contactSet_;
  mc_rtc::duration_ms dt_compute_minkSum_;
  mc_rtc::duration_ms dt_update_planes_;
  mc_rtc::duration_ms dt_compute_guiTrianglesContacts_;
  mc_rtc::duration_ms dt_compute_guiTrianglesRegions_ = mc_rtc::duration_ms::zero();
  mc_rtc::duration_ms dt_zeroMoment_intersection_;
  std::map<std::string, ContactTimers> contactsTimers_;

  // map of polytope triangles for display
  std::map<std::string, std::vector<std::array<Eigen::Vector3d, 3>>> frictionConesTrianglesMap_;
  std::map<std::string, std::vector<std::array<Eigen::Vector3d, 3>>> forcePolyTrianglesMap_;
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
