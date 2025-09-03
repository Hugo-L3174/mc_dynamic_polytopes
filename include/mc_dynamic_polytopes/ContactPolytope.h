#pragma once
#include <mc_dynamic_polytopes/AsyncJob.h>

#include <chrono>
#include <future>

#include <SpaceVecAlg/SpaceVecAlg>
#include <boost/shared_ptr.hpp>
#include <boost/smart_ptr/make_shared_object.hpp>

#include <mc_rtc/clock.h>
#include <mc_rtc/gui.h>
#include <mc_rtc/logging.h>

#include <mc_dynamic_polytopes/GUIComputations.h>
#include <mc_dynamic_polytopes/PolytopeFuncs.h>
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
constexpr double defaultForceScale = 1;

struct ContactTimers
{
  std::atomic<double> dt_frictionCone = 0.;
  std::atomic<double> dt_forcePolytope = 0.;
  std::atomic<double> dt_intersection = 0.;
  std::atomic<double> dt_contactTotal = 0.;
  std::atomic<double> dt_copyInputs = 0.;
  std::atomic<double> dt_constructor = 0.;
  std::atomic<double> dt_double_description = 0.;
  std::atomic<double> dt_qhull = 0.;
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

  // number of sides for cones linearization
  std::atomic<int> numberOfFrictionSides = defaultFrictionSides;
  // force scaling factor (alpha to be used to transfer between contacts)
  std::atomic<double> forceScalingFactor = defaultForceScale;
  // friction coeff for the cones
  double frictionCoefficient = defaultFrictionCoeff;
};

struct ContactPolytopeJob : public MakeAsyncJob<ContactPolytopeJob, ContactPolytopeInput, ContactPolytopeResult>
{
  // Options for easier display of the computation steps
  bool combineWithFriction_ = true;
  bool DDfrictionCones_ = false;
  bool HrepMode_ = false;
  double guiScale_ = 0.001;

  // Polyhedron configs
  mc_rtc::gui::PolyhedronConfig polyForceConfig_;
  mc_rtc::gui::PolyhedronConfig polyMomentConfig_;

  // Additional timers for polytope computation steps
  ContactTimers contactTimers;

  mc_rbdyn::Contact * contactRBDyn_; // for gui friction coeff

  // Actual computation
  void buildActuationPolytopeFromContact(const ContactPolytopeInput & input,
                                         boost::shared_ptr<Polytope_Rn> & actuationPolytope,
                                         double forceScalingFactor,
                                         unsigned int dim);

  // CRTP actual computation
  ContactPolytopeResult computeJob();

  // CRTP: deferred logger implementation
  void addToLoggerImpl()
  {
    auto contact = input_.contactName;
    auto contactPrefix = "perf_" + loggerPrefix_ + "_" + contact + "_async_";
    logger_->addLogEntry(contactPrefix + "frictionCone [ms]", this,
                         [this]() { return contactTimers.dt_frictionCone.load(); });
    logger_->addLogEntry(contactPrefix + "forcePolytope [ms]", this,
                         [this]() { return contactTimers.dt_forcePolytope.load(); });
    logger_->addLogEntry(contactPrefix + "intersection [ms]", this,
                         [this]() { return contactTimers.dt_intersection.load(); });
    logger_->addLogEntry(contactPrefix + "contactTotal [ms]", this,
                         [this]() { return contactTimers.dt_contactTotal.load(); });
  }

  // CRTP: deferred GUI implementation
  void addToGUIImpl()
  {
    const auto & result = *lastResult_;
    auto coeffsCat = guiCategory_;
    auto contactsCat = guiCategory_;
    contactsCat.push_back("Contact Polytopes");
    auto contact = input_.contactName;
    using namespace mc_rtc::gui;

    gui_->addElement(this, contactsCat,
                     Polyhedron(fmt::format(contact + " frictions"), polyForceConfig_,
                                [&result]() { return result.frictionConeTriangles; }),
                     Polyhedron(fmt::format(contact + " forces"), polyMomentConfig_,
                                [&result]() { return result.forcePolyTriangles; }),
                     Polyhedron(fmt::format(contact + " moments"), polyMomentConfig_,
                                [&result]() { return result.momentPolytopesTriangles; }));

    gui_->addElement(this, coeffsCat,
                     NumberSlider(
                         fmt::format(contact + " force alpha [0.001-1]"),
                         [this, contact]() { return input_.forceScalingFactor.load(); },
                         [this](double scale) { input_.forceScalingFactor = scale; }, 0.001, 1.0),
                     IntegerInput(
                         fmt::format(contact + " number of friction sides"),
                         [this]() { return input_.numberOfFrictionSides.load(); },
                         [this](int nbFrictionSides) { input_.numberOfFrictionSides = nbFrictionSides; }),
                     NumberInput(
                         fmt::format(contact + " friction coefficient"), [this]() { return contactRBDyn_->friction(); },
                         [this](double mu) { contactRBDyn_->friction(mu); }));
  }

  void load(const mc_rtc::Configuration & config, const std::string & contactName)
  {
    auto loadPolyConfigs = [&](const mc_rtc::Configuration & conf)
    {
      conf("polyhedronForce", polyForceConfig_);
      conf("polyhedronMoment", polyMomentConfig_);
    };

    loadPolyConfigs(config);
    if(auto contactsConf = config.find("contacts"))
    {
      if(auto contactConf = contactsConf->find(contactName))
      {
        loadPolyConfigs(*contactConf);
      }
    }
  }
};

} // namespace mc_dynamic_polytopes
