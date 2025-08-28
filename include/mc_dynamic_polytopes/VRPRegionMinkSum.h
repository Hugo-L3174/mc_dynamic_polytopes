#pragma once

#include <mc_rbdyn/Robot.h>
#include <atomic>
#include <mc_dynamic_polytopes/AsyncJob.h>
#include <mc_dynamic_polytopes/ContactPolytope.h>
#include <politopix/politopixAPI.h>

namespace mc_dynamic_polytopes
{

struct VRPMinkSumTimers
{
  std::atomic<double> dt_update_planes{0};
  std::atomic<double> dt_zeroMoment_intersection{0};
  std::atomic<double> dt_compute_minkSum{0};
  std::atomic<double> dt_compute_guiTrianglesRegions_{0};
};

struct VRPRegionMinkSumJobInput
{
  /**
   * Partial initialization from the robot instance
   * You need to instanciate the remaining members yourself
   */
  void initialize_robot(const mc_rbdyn::Robot & robot)
  {
    robotMass = robot.mass();
    comPosition = robot.com();
  }

  boost::shared_ptr<Polytope_Rn> zmpRegion;
  std::map<std::string, ContactPolytopeResult> contactsPolytopes;
  std::map<std::string, sva::PTransformd> contactsPose;
  double robotMass = 0.0;
  Eigen::Vector3d comPosition = Eigen::Vector3d::Zero();
};

struct VRPRegionMinkSumJobResult
{
  boost::shared_ptr<Polytope_Rn> CWCForces;
  boost::shared_ptr<Polytope_Rn> CWCMoments;
  boost::shared_ptr<Polytope_Rn> zeroMomentRegion;
  // Internal matrices of planes and offsets of the regions for constraints
  std::pair<Eigen::MatrixXd, Eigen::VectorXd> DCMVRPPlanes; // Matrix constraint for force polytope
                                                            //
                                                            //
  // internal normals matrix and offsets vector of the zero moment region (subset of the DCM region) for QP constraint
  // or check
  std::pair<Eigen::MatrixXd, Eigen::VectorXd> zeroMomentPlanes; // Matrix constraint for zero moment region
  std::vector<std::array<Eigen::Vector3d, 3>> CWCForceTriangles;
  std::vector<std::array<Eigen::Vector3d, 3>> CWCMomentTriangles;
  std::vector<std::array<Eigen::Vector3d, 3>> ZMPTriangles;
  std::vector<std::array<Eigen::Vector3d, 3>> zeroMomentTriangles;
};

struct VRPRegionMinkSumJob : MakeAsyncJob<VRPRegionMinkSumJob, VRPRegionMinkSumJobInput, VRPRegionMinkSumJobResult>
{
  bool withVRPOffset_ = true;
  bool withMoments_ = false;
  double guiScale_ = 0.001;

protected:
  VRPMinkSumTimers vrpTimers_;
  VRPRegionMinkSumJobResult result_;
  std::string name_ = "VRPRegionMinkSumJob";
  mc_rtc::gui::PolyhedronConfig polyForceConfig_;
  mc_rtc::gui::PolyhedronConfig polyMomentConfig_;
  mc_rtc::gui::PolyhedronConfig polyZMPConfig_;
  mc_rtc::gui::PolyhedronConfig polyZeroMomentAreaConfig_;

public: // XXX: should this be public?
  /**
   * Computes the mink sum of the contact regions, scale to get eCMP region, intersect with ZMP region, and translate
   * everything to get VRP region */
  VRPRegionMinkSumJobResult computeJob(); // void computeVRPRegionWithMinkSum();

  void addToLoggerImpl()
  {
    auto prefix = "perf_" + loggerPrefix_ + "_async_";
    logger_->addLogEntry(prefix + "compute_minkSum [ms]", this,
                         [this]() { return vrpTimers_.dt_compute_minkSum.load(); });
    logger_->addLogEntry(prefix + "zeroMoment_intersection [ms]", this,
                         [this]() { return vrpTimers_.dt_zeroMoment_intersection.load(); });
    logger_->addLogEntry(prefix + "update_planes [ms]", this, [this]() { return vrpTimers_.dt_update_planes.load(); });
    logger_->addLogEntry(prefix + "compute_guiTrianglesRegions [ms]", this,
                         [this]() { return vrpTimers_.dt_compute_guiTrianglesRegions_.load(); });
  }

  void addToGUIImpl()
  {
    auto CWCCat = guiCategory_;
    CWCCat.push_back("Contact Wrench Cone");

    gui_->addElement(
        this, CWCCat,
        mc_rtc::gui::Polyhedron("CWC forces", polyForceConfig_, [this]() { return result_.CWCForceTriangles; }),
        mc_rtc::gui::Polyhedron("CWC moments", polyMomentConfig_, [this]() { return result_.CWCMomentTriangles; }),
        mc_rtc::gui::Polyhedron("ZMP area", polyZMPConfig_, [this]() { return result_.ZMPTriangles; }),
        mc_rtc::gui::Polyhedron("Zero moment region", polyZeroMomentAreaConfig_,
                                [this]() { return result_.zeroMomentTriangles; }));
  }

  void load(const mc_rtc::Configuration & config)
  {
    config("polyhedronForce", polyForceConfig_);
    config("polyhedronMoment", polyMomentConfig_);
    config("polyhedronZMP", polyZMPConfig_);
    config("polyhedronZeroMomentArea", polyZeroMomentAreaConfig_);
  }

protected:
  // Computes the minkowsky sum of the given friction cones and puts the result in the CWC_ polytope object
  void computeMinkowskySumPolitopix();
  // Computes the eCMP region from the CWC
  // /!\ this modifies the CWC to invert and scale it, the result is still in the CWC polytope !
  void computeECMPRegion(const Eigen::Vector3d & comPosition, double robotMass);
  /* Computes the intersection between the eCMP region and the ZMP region to get the zero moment region,
  put in zeroMomentRegion_
  */
  void computeZeroMomentIntersection();

  // Translates the eCMP region and the zero-moment region with a vertical Delta z offset to get actual DCM region
  void VRPtranslation(double deltaZ);

  // Updates the internal maps of triangles of the total regions for gui display
  void updateTrianglesRegionsGUIPolitopix();
};

} // namespace mc_dynamic_polytopes
