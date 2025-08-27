#pragma once

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <politopix/Voronoi_Rn.h>
#include <thread>

#include <RBDyn/Coriolis.h>
#include <SpaceVecAlg/SpaceVecAlg>
#include <boost/smart_ptr/make_shared_object.hpp>

#include <mc_control/fsm/Controller.h>
#include <mc_rtc/clock.h>
#include <mc_rtc/gui.h>
#include <mc_rtc/logging.h>
#include <Tasks/QPContacts.h>

#include <mc_dynamic_polytopes/ContactPolytope.h>
#include <mc_dynamic_polytopes/ZMPRegion.h>
#include <politopix/PolyhedralAlgorithms_Rn.h>
#include <politopix/PrismaticPolyhedron_Rn.h>
#include <politopix/politopixAPI.h>

#include <libqhull_r/qhull_ra.h>

namespace mc_dynamic_polytopes
{

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
    if(minkSumThread_.joinable())
    {
      minkSumThread_.join();
    }
  }

  // Set current contact set to be used for next computation: contact name and contact reference pose, ie the frame of
  // the target surface of the contact (only the orientation is used for the friction cone)
  void setControllerContacts(const std::map<std::string, sva::PTransformd> & contacts)
  {
    // XXX remove this map from header if switch to rbdyn contacts
    refContactTransforms_.clear();
    for(const auto & contact : contacts)
    {
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

  void setControllerContacts(const std::map<std::string, mc_rbdyn::Contact &> & contacts)
  {
    refContactTransforms_.clear();
    contactsRBDyn_ = contacts;
    for(const auto & contact : contacts)
    {
      const auto & contactName = contact.first;
      const auto & contactObj = contact.second;
      // XXX X_r_r is not updated automatically by default, it needs to be done in the controller
      // FIXME identity matrix is an assumption that the contact plane is the controlled frame for friction
      if(contactName == contactObj.r1Surface()->name())
      {
        refContactTransforms_.emplace(contactName, sva::PTransformd::Identity()); // contactObj.X_r2s_r1s().inv());
      }
      else
      {
        // if second robot, assume it is targeted and put identity
        refContactTransforms_.emplace(contactName, sva::PTransformd::Identity());
      }
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

  void addToGUI(mc_rtc::gui::StateBuilder * gui);
  void addToLogger(mc_rtc::Logger & logger);
  void removeFromLogger(mc_rtc::Logger & logger);

  // ------------------------------------------------------> public computation functions

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

  // ------------------------------------------------------> Plane constraints getters

  // Returns the internal normals matrix and offsets vector of the eCMP region for QP constraint or check
  const HRepXd & getVRPPlanes()
  {
    std::lock_guard<std::mutex> lock(VRPPlanesMutex_);
    return DCMVRPPlanes_;
  };

  // Returns the internal normals matrix and offsets vector of the zero moment region (subset of the DCM region) for QP
  // constraint or check
  const HRepXd & getZeroMomentPlanes()
  {
    std::lock_guard<std::mutex> lock(zeroMomentPlanesMutex_);
    return zeroMomentPlanes_;
  };

  // Returns the desired cone H representation
  // Prefer using the force polytope planes if they are computed
  const HRepXd & getFrictionConePlanes(const std::string & contactName)
  {
    if(feasiblePolytopesJobs_.count(contactName))
    {
      auto & job = feasiblePolytopesJobs_[contactName];
      if(auto result = job.lastResult())
      {
        return result->frictionConesPlanes;
      }
    }
    mc_rtc::log::error_and_throw("Cannot get friction cone planes for contact {}", contactName);
  };

  // Returns the desired contact force polytope H representation
  const HRepXd & getForcePolyPlanes(const std::string & contactName)
  {
    if(feasiblePolytopesJobs_.count(contactName))
    {
      auto & job = feasiblePolytopesJobs_[contactName];
      if(auto result = job.lastResult())
      {
        return result->forcePolyPlanes;
      }
    }
    mc_rtc::log::error_and_throw("Cannot get force poly planes for contact {}", contactName);
  };

  const auto & robot() const noexcept
  {
    return robot_;
  }

  /**
   * XXX: MUST be called every iteration
   */
  void computeRegions();

  // Projects the given point on the VRP region. Returns the given point if it is already inside
  // Eigen::Vector3d projectPointInVRPRegion(Eigen::Vector3d testedPoint);
  //
  // Eigen::Vector3d projectPointInZeroMomentRegion(Eigen::Vector3d testedPoint);

protected:
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
  boost::shared_ptr<Polytope_Rn> computeZMPRegion(Eigen::Vector3d comPosition);

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

  // Updates the feasible polytope representation internal to rbdyn contacts
  void updateRBDynPolytopes(const Eigen::MatrixXd & Normals,
                            const Eigen::VectorXd & Offsets,
                            mc_rbdyn::Contact & contactRBDyn);

  // From the current contact set, deduce what contacts need to be removed from computation compared to last iteration
  void setCurrentContacts()
  {
    // take the previous set of contacts
    contactsToRemove_ = activeContacts_;
    auto previousActiveContacts = activeContacts_;
    activeContacts_.clear();
    for(const auto & contact : refContactTransforms_)
    {
      // add every contact given to the active contacts
      activeContacts_.emplace(contact.first);
      // every active contact does not need to be removed
      contactsToRemove_.erase(contact.first);
    }
  }

  // Eigen::Vector3d projectPointInPolytope(Eigen::Vector3d testedPoint, boost::shared_ptr<Polytope_Rn> & polytope);

  // sanity check
  bool checkGravityCenterInPolytope(boost::shared_ptr<Polytope_Rn> & polytope);

  // ------------------------------------------------------> GUI getters

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

  // ------------------------------------------------------> Internal variables

  std::string name_;
  std::vector<std::string> guiCategory_ = {"DynamicPolytopes"};
  const mc_rbdyn::Robot & robot_;
  mc_rtc::gui::StateBuilder * gui_ = nullptr;
  mc_rtc::Logger * logger_ = nullptr;
  mc_rtc::Configuration config_;
  std::set<std::string> possibleContacts_;
  std::set<std::string> activeContacts_;
  bool contactsUpdated_ = false;
  std::set<std::string> contactsToRemove_;

  // map of the contact transforms X_r1_r2, with r1 the controlled robot frame and r2 the target frame of the contact
  // This is used to compute the orientation of the friction cone as it should be in the r1 frame (same as the force
  // polytope) for distribution, and not in world frame
  std::map<std::string, sva::PTransformd> refContactTransforms_;

  std::map<std::string, mc_rbdyn::Contact &> contactsRBDyn_;

  bool withMoments_;
  bool computeRegions_;
  bool HrepMode_;
  bool withVRPOffset_ = true;
  bool combineWithFriction_ = true;
  bool DDfrictionCones_ = false;

  // politopix

  std::map<std::string, ContactPolytopeJob>
      feasiblePolytopesJobs_; /// For each contact store a job to compute its feasibility polytope asynchronously

  boost::shared_ptr<Polytope_Rn> CWCForces_;
  boost::shared_ptr<Polytope_Rn> CWCMoments_;
  // boost::shared_ptr<Polytope_Rn> eCMPRegion_;

  boost::shared_ptr<Polytope_Rn> zmpRegion_;
  boost::shared_ptr<Polytope_Rn> zeroMomentRegion_;

  // threading things
  // atomic bool as condition to keep computing in loop
  std::atomic<bool> computing_{false};
  // condition variable to signal start of loop for main thread
  std::condition_variable cv_;
  std::mutex contactSetMutex_;
  ZMPRegionJob zmpRegionJob_;
  std::thread minkSumThread_;
  std::atomic<bool> VRPRegionComputed_{true};
  std::mutex VRPPlanesMutex_;
  std::mutex zeroMomentPlanesMutex_;

  std::mutex CWCMutex_;
  std::mutex ZMPMutex_;
  std::mutex zeroMomentMutex_;

  // gui mutexes
  std::mutex CWCForceTrianglesMutex_;
  std::mutex CWCMomentTrianglesMutex_;
  std::mutex ZMPTrianglesMutex_;
  std::mutex zeroMomentTrianglesMutex_;

  // Internal matrices of planes and offsets of the regions for constraints
  HRepXd DCMVRPPlanes_;

  HRepXd zeroMomentPlanes_;

  // timers to measure computation times
  mc_rtc::duration_ms dt_loop_total_;
  mc_rtc::duration_ms dt_compute_contactSet_;
  mc_rtc::duration_ms dt_compute_minkSum_;
  mc_rtc::duration_ms dt_update_planes_;
  mc_rtc::duration_ms dt_compute_guiTrianglesContacts_;
  mc_rtc::duration_ms dt_compute_guiTrianglesRegions_ = mc_rtc::duration_ms::zero();
  mc_rtc::duration_ms dt_zeroMoment_intersection_;
  mc_rtc::duration_ms dt_debug_;

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

} // namespace mc_dynamic_polytopes
