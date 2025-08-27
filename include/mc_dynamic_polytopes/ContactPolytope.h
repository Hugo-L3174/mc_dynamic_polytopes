#pragma once

#include <chrono>
#include <future>

#include <SpaceVecAlg/SpaceVecAlg>
#include <boost/shared_ptr.hpp>
#include <boost/smart_ptr/make_shared_object.hpp>

#include <mc_rtc/clock.h>
#include <mc_rtc/gui.h>
#include <mc_rtc/logging.h>

#include <mc_dynamic_polytopes/GUIComputations.h>
#include <mc_dynamic_polytopes/WrenchCones.h>

#include <politopix/PolyhedralAlgorithms_Rn.h>
#include <politopix/PrismaticPolyhedron_Rn.h>
#include <politopix/politopixAPI.h>

#include <libqhull_r/qhull_ra.h>

namespace mc_dynamic_polytopes
{

using HRepXd = std::pair<Eigen::MatrixXd, Eigen::VectorXd>;
constexpr int defaultFrictionSides = 5;
constexpr double defaultFrictionCoeff = 0.5;

struct ContactTimers
{
  using duration_ms = mc_rtc::duration_ms;
  duration_ms dt_frictionCone = duration_ms::zero();
  duration_ms dt_forcePolytope = duration_ms::zero();
  duration_ms dt_intersection = duration_ms::zero();
  duration_ms dt_contactTotal = duration_ms::zero();
  duration_ms dt_startAsync = duration_ms::zero();
  duration_ms dt_copyInputs = duration_ms::zero();
  duration_ms dt_constructor = duration_ms::zero();
  duration_ms dt_total_job = duration_ms::zero();
};

struct ContactPolytopeResult
{
  ContactPolytopeResult()
  {
    frictionCone = boost::make_shared<Polytope_Rn>();
    actuationPolytope = boost::make_shared<Polytope_Rn>();
    frictionConeMoments = boost::make_shared<Polytope_Rn>();
    // XXX: what's the point of this initialization
    // Here initialize planes as simple friction cones to have a sane constraint at first iteration
    // The rotX_r1_r2 matrix would be the identity in a default case
    HRepXd newPlanes;
    newPlanes.first =
        generatePolyhedralConeHRep(defaultFrictionSides, Eigen::Matrix3d::Identity(), defaultFrictionCoeff);
    newPlanes.second = Eigen::VectorXd::Zero(newPlanes.first.rows());
    frictionConesPlanes = newPlanes;
    forcePolyPlanes = newPlanes;
  }

  void updateTrianglesGUIPolytopix(double guiScale, const sva::PTransformd & contactPose)
  {
    if(actuationPolytope->dimension() == 3)
    {
      update3DPolyTrianglesPolitopix(frictionCone, frictionConeTriangles, guiScale, contactPose);
      update3DPolyTrianglesPolitopix(actuationPolytope, forcePolyTriangles, guiScale, contactPose);
    }
    else
    {
      update6DPolyTrianglesPolitopix(actuationPolytope, forcePolyTriangles, momentPolytopesTriangles, guiScale,
                                     contactPose);
    }
  }

  std::string contactName;
  // We keep the friction cones as polytope objects, but they are actually polyhedral cones and will not be bounded
  boost::shared_ptr<Polytope_Rn> frictionCone;
  // The bounded actuation polytopes
  boost::shared_ptr<Polytope_Rn> actuationPolytope;
  boost::shared_ptr<Polytope_Rn> frictionConeMoments;

  HRepXd frictionConesPlanes;
  HRepXd forcePolyPlanes;

  std::vector<std::array<Eigen::Vector3d, 3>> frictionConeTriangles;
  std::vector<std::array<Eigen::Vector3d, 3>> forcePolyTriangles;
  std::vector<std::array<Eigen::Vector3d, 3>> momentPolytopesTriangles;
};

struct ContactPolytopeInput
{
  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  std::shared_ptr<mc_rbdyn::Surface> surface;
  sva::MotionVecd accW;
  std::vector<std::vector<double>> tl, tu;
  sva::PTransformd surfacePose;

  // TODO: sane defaults
  std::string contactName;
  sva::PTransformd refContactTransform;
  int numberOfFrictionSides;
  double forceScalingFactors;
  double frictionCoefficients;
};

struct ContactPolytopeJob
{
  ~ContactPolytopeJob()
  {
    removeFromLogger();
    removeFromGUI();
  }

  ContactPolytopeInput input;
  ContactTimers timers;

  // XXX: what to do about these options?
  // Options for easier display of the computation steps
  // TODO: turn into a schema
  bool combineWithFriction_ = true;
  bool DDfrictionCones_ = false;
  bool HrepMode_ = false;
  double guiScale_ = 0.001;

  /**
   * Starts the async job defined in \ref computePolytopeJob()
   * One must set the input before this
   */
  void startAsync()
  {
    auto start_async = mc_rtc::clock::now();
    running_ = true;
    futurePolytope = std::async(std::bind(&ContactPolytopeJob::computePolytopeJob, this));
    timers.dt_startAsync = mc_rtc::clock::now() - start_async;
  }

  /**
   * Checks whether the async task has been started
   */
  bool running() const noexcept
  {
    return running_;
  }

  /**
   * Checks whether the job has completed and stores it to lastResult_
   *
   * \returns true when the result has been computed, false otherwise
   */
  bool checkResult()
  {
    if(futurePolytope.wait_for(std::chrono::seconds(0)) == std::future_status::ready)
    {
      lastResult_ = futurePolytope.get();
      if(logger_ && !inLogger_)
      {
        addToLogger_();
        inLogger_ = true;
      }
      if(gui_ && !inGUI_)
      {
        addToGUI_();
        inGUI_ = true;
      }
      running_ = false;
      return true;
    }
    return false;
  }

  /**
   * Get the last result
   * One should check if we have a result using checkResult
   */
  const std::optional<ContactPolytopeResult> & lastResult() const noexcept
  {
    return lastResult_;
  }

  /**
   * Add job related information to the logger
   * This is defered until:
   * - addToLogger has been called
   * - the first async job has been completed
   */
  void addToLogger(mc_rtc::Logger & logger, const std::string & prefix)
  {
    if(logger_) return;
    logger_ = &logger;
    auto contact = input.contactName;
    auto contactPrefix = "perf_" + prefix + "_" + contact + "_";
    addToLogger_ = [this, &logger, contactPrefix]()
    {
      logger.addLogEntry(contactPrefix + "frictionCone", this, [this]() { return timers.dt_frictionCone.count(); });
      logger.addLogEntry(contactPrefix + "forcePolytope", this, [this]() { return timers.dt_forcePolytope.count(); });
      logger.addLogEntry(contactPrefix + "intersection", this, [this]() { return timers.dt_intersection.count(); });
      logger.addLogEntry(contactPrefix + "Total", this, [this]() { return timers.dt_contactTotal.count(); });
      logger.addLogEntry(contactPrefix + "startAsync", this, [this]() { return timers.dt_startAsync.count(); });
      logger.addLogEntry(contactPrefix + "copyInputs", this, [this]() { return timers.dt_copyInputs.count(); });
      logger.addLogEntry(contactPrefix + "constructor", this, [this]() { return timers.dt_constructor.count(); });
      logger.addLogEntry(contactPrefix + "totalJob", this, [this]() { return timers.dt_total_job.count(); });
    };
  }

  /**
   * Add job related information to the GUI
   * This is defered until:
   * - addToGUI has been called
   * - the first async job has been completed
   */
  void addToGUI(mc_rtc::gui::StateBuilder & gui, const std::vector<std::string> & category)
  {
    if(gui_) return;
    gui_ = &gui;
    addToGUI_ = [this, &gui, category]()
    {
      const auto & result = *lastResult_;
      auto coeffsCat = category;
      auto contactsCat = category;
      contactsCat.push_back("Contact Polytopes");
      auto contact = input.contactName;

      gui.addElement(this, contactsCat,
                     mc_rtc::gui::Polyhedron(fmt::format(contact + " frictions"), polyForceConfig_,
                                             [&result]() { return result.frictionConeTriangles; }),
                     mc_rtc::gui::Polyhedron(fmt::format(contact + " forces"), polyMomentConfig_,
                                             [&result]() { return result.forcePolyTriangles; }),
                     mc_rtc::gui::Polyhedron(fmt::format(contact + " moments"), polyMomentConfig_,
                                             [&result]() { return result.momentPolytopesTriangles; }));
    };
  }

  void removeFromLogger()
  {
    if(logger_)
    {
      logger_->removeLogEntries(this);
    }
  }

  void removeFromGUI()
  {
    if(gui_)
    {
      gui_->removeElements(this);
    }
  }

  void load(const mc_rtc::Configuration & config)
  {
    // XXX: simplify with schema
    config("gui")("polyhedronForce")("triangle_color", polyForceConfig_.triangle_color);
    config("gui")("polyhedronForce")("show_triangle", polyForceConfig_.show_triangle);
    config("gui")("polyhedronForce")("use_triangle_color", polyForceConfig_.use_triangle_color);
    config("gui")("polyhedronForce")("edges", polyForceConfig_.edge_config);
    config("gui")("polyhedronForce")("show_edges", polyForceConfig_.show_edges);
    config("gui")("polyhedronForce")("fixed_edge_color", polyForceConfig_.fixed_edge_color);
    config("gui")("polyhedronForce")("vertices")("color", polyForceConfig_.vertices_config.color);
    config("gui")("polyhedronForce")("vertices")("scale", polyForceConfig_.vertices_config.scale);
    config("gui")("polyhedronForce")("show_vertices", polyForceConfig_.show_vertices);
    config("gui")("polyhedronForce")("fixed_vertices_color", polyForceConfig_.fixed_vertices_color);

    config("gui")("polyhedronMoment")("triangle_color", polyMomentConfig_.triangle_color);
    config("gui")("polyhedronMoment")("show_triangle", polyMomentConfig_.show_triangle);
    config("gui")("polyhedronMoment")("use_triangle_color", polyMomentConfig_.use_triangle_color);
    config("gui")("polyhedronMoment")("edges", polyMomentConfig_.edge_config);
    config("gui")("polyhedronMoment")("show_edges", polyMomentConfig_.show_edges);
    config("gui")("polyhedronMoment")("fixed_edge_color", polyMomentConfig_.fixed_edge_color);
    config("gui")("polyhedronMoment")("vertices")("color", polyMomentConfig_.vertices_config.color);
    config("gui")("polyhedronMoment")("vertices")("scale", polyMomentConfig_.vertices_config.scale);
    config("gui")("polyhedronMoment")("show_vertices", polyMomentConfig_.show_vertices);
    config("gui")("polyhedronMoment")("fixed_vertices_color", polyMomentConfig_.fixed_vertices_color);
  }

protected: // bookeeping for the async job
  bool running_ = false;
  bool inLogger_ = false;
  bool inGUI_ = false;
  std::future<ContactPolytopeResult> futurePolytope;
  std::optional<ContactPolytopeResult> lastResult_;
  ContactPolytopeResult computePolytopeJob();

protected: // deferred loggin/GUI calls
  mc_rtc::Logger * logger_ = nullptr;
  std::function<void()>
      addToLogger_; ///< Actual logging implementation to be called after the first result is available

  mc_rtc::gui::StateBuilder * gui_ = nullptr;
  std::function<void()> addToGUI_; ///< Actual gui implementation to be called after the first result is available
  mc_rtc::gui::PolyhedronConfig polyForceConfig_;
  mc_rtc::gui::PolyhedronConfig polyMomentConfig_;

  // actual computation
protected:
  // compute the force polytope of the contact from the wrench limits of the limb actuating it
  // now also taking a scale factor (0-1) for force
  void buildActuationPolytopeFromContact(const ContactPolytopeInput & input,
                                         boost::shared_ptr<Polytope_Rn> & actuationPolytope,
                                         double forceScalingFactor,
                                         unsigned int dim);
};

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
