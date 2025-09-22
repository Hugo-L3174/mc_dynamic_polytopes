#include "mc_dynamic_polytopes/Time.h"
#include <mc_dynamic_polytopes/DynamicPolytope.h>

namespace mc_dynamic_polytopes
{

DynamicPolytope::DynamicPolytope(const std::string & name,
                                 const mc_rbdyn::Robot & robot,
                                 const mc_rtc::Configuration & dynamicPolyConfig)
: name_(fmt::format("DynamicPolytope_" + name)), guiCategory_({"DynamicPolytopes", name}), robot_(robot),
  config_(dynamicPolyConfig)
{
  withMoments_ = dynamicPolyConfig("withMoments", false);
  computeRegions_ = dynamicPolyConfig("computeRegions", true);
  HrepMode_ = dynamicPolyConfig("HrepMode", false);
  combineWithFriction_ = dynamicPolyConfig("withFriction", true);

  // Init dimension
  unsigned int dim = withMoments_ ? 6 : 3;
  Rn::setDimension(dim);
  Rn::setTolerance(1e-07);
}

// Check contacts from ctl
//
// Start async thread for each new contact
//    For existing threads, check if computation is finished (boolean)
//    lock inputstate
//    copy inputstate
//    unlock inputstate
//
//    work on copy of inputstate
//    ... thread computes ...
//    setcomputation finished flag
//
// Stop threads for removed contacts
//
// Get result from all threads
// Update contact in ctl
void DynamicPolytope::computeRegions()
{
  auto start_loop = mc_rtc::clock::now();
  withMoments_ ? Rn::setDimension(6) : Rn::setDimension(3);

  // Get the set of active contacts and contacts to remove
  setCurrentContacts();

  auto & robot = robot_;
  // Step 1.1: launch individual force polytopes and friction cones calculations in separate threads as they are
  // independant, then their intersection
  for(const auto & contactName : activeContacts_)
  {
    // TODO: Build input for each contact
    // and replace this call
    auto start_constructor = mc_rtc::clock::now();
    auto & job = feasiblePolytopesJobs_[contactName];
    job.contactRBDyn_ = &contactsRBDyn_.at(contactName);
    if(!job.startedOnce())
    {
      job.load(config_("gui", mc_rtc::Configuration{}),
               contactName); // XXX should load per-contact config instead of global one
    }
    job.contactTimers.dt_constructor = mc_rtc::elapsed_ms_count(start_constructor);
    job.HrepMode_ = HrepMode_; // XXX

    if(!job.running())
    {
      // Initialize inputs required by the contact polytope computation
      // It is only safe to modify job.input() while the async task is not running
      job.contactTimers.dt_copyInputs = mc_rtc::timedExecution(
                                            [&]()
                                            {
                                              auto & input = job.input();
                                              input.contactName = contactName;
                                              input.mb = robot.mb();
                                              input.mbc = robot.mbc();
                                              input.tl = robot.tl();
                                              input.tu = robot.tu();
                                              input.surface = robot.surface(contactName).copy();
                                              input.surfacePose = input.surface->X_0_s(robot);
                                              input.accW = robot.accW();
                                              input.refContactTransform = refContactTransforms_.at(contactName);
                                              input.frictionCoefficient = contactsRBDyn_.at(contactName).friction();
                                            })
                                            .count();

      // Job starts as an async task, use job.checkResult() later to know whether it is finished and retrive its value
      job.startAsync();
      if(logger_)
      {
        job.addToLogger(*logger_, name_ + "_ContactPolytope_" + contactName);
      }
      if(gui_)
      {
        job.addToGUI(*gui_, guiCategory_);
      }
    }
  }

  // Step 2: wait for finished individual feasible regions
  // Wait for finished threads and join them, then update their matrix constraints
  dt_update_rbdyn_ = mc_rtc::timedExecution(
      [&, this]()
      {
        for(const auto & contactName : activeContacts_)
        {
          auto & job = feasiblePolytopesJobs_[contactName];
          if(job.checkResult())
          {
            const auto & result = *job.lastResult();
            updateRBDynPolytopes(result.forcePolyPlanes.first, result.forcePolyPlanes.second,
                                 contactsRBDyn_.at(contactName));
          }
        }
      });

  // Remove contacts
  for(const auto & removeContact : contactsToRemove_)
  {
    feasiblePolytopesJobs_.erase(removeContact);
  }

  // Step 1.2: launch ZMP region calculation at the same time, also independant
  if(computeRegions_)
  {
    if(!zmpRegionJob_.running())
    {
      zmpRegionJob_.input().initialize(robot, activeContacts_);
      zmpRegionJob_.startAsync();

      if(logger_)
      {
        zmpRegionJob_.addToLogger(*logger_, name_ + "_ZMPRegion");
      }
      if(gui_)
      {
        auto zmpRegionCategory = guiCategory_;
        zmpRegionCategory.push_back("ZMP Region");
        zmpRegionJob_.addToGUI(*gui_, zmpRegionCategory);
      }
    }
    zmpRegionJob_.checkResult();
  }

  dt_compute_contactSet_ = mc_rtc::elapsed_ms(start_loop);

  // Steps 3-6: launch the rest everytime the previous full region was computed
  if(computeRegions_)
  {
    // After checking results for feasiblePolytopesJobs_ and zmpRegionJob_
    bool allContactsReady = std::all_of(activeContacts_.begin(), activeContacts_.end(),
                                        [&](const auto & contactName)
                                        {
                                          return feasiblePolytopesJobs_.count(contactName)
                                                 && feasiblePolytopesJobs_.at(contactName).lastResult().has_value();
                                        });
    bool zmpReady = zmpRegionJob_.lastResult().has_value();

    if(allContactsReady && zmpReady && !VRPRegionMinkSumJob_.running())
    {
      if(!VRPRegionMinkSumJob_.startedOnce())
      {
        VRPRegionMinkSumJob_.load(config_("gui", mc_rtc::Configuration{}));
      }
      auto & input = VRPRegionMinkSumJob_.input();
      input.initialize_robot(robot_);
      input.contactsPolytopes.clear();
      input.contactsPose.clear();

      for(const auto & contactName : activeContacts_)
      {
        const auto & job = feasiblePolytopesJobs_.at(contactName);
        input.contactsPolytopes[contactName] = *job.lastResult();
        input.contactsPose[contactName] = refContactTransforms_.at(contactName);
      }

      input.zmpRegion = zmpRegionJob_.lastResult()->zmpRegion;

      VRPRegionMinkSumJob_.startAsync();
      if(logger_)
      {
        VRPRegionMinkSumJob_.addToLogger(*logger_, name_ + "_VRPRegionMinkSum");
      }
      if(gui_)
      {
        auto vrpRegionCategory = guiCategory_;
        vrpRegionCategory.push_back("VRP Region");
        VRPRegionMinkSumJob_.addToGUI(*gui_, vrpRegionCategory);
      }
    }

    if(VRPRegionMinkSumJob_.checkResult())
    {
      // const auto & result = *VRPRegionMinkSumJob_.lastResult();
      // Use result.CWCForces, result.CWCMoments, result.zeroMomentRegion, etc.
    }
  }

  // Regions GUI is now updated at the end of their thread to ensure scaling and translation are done before updati
  dt_loop_total_ = mc_rtc::elapsed_ms(start_loop);
}

void DynamicPolytope::buildFrictionConeFromContactWithVrep(int numberOfFrictionSides,
                                                           const sva::PTransformd contactSurface,
                                                           boost::shared_ptr<Polytope_Rn> & frictionCone,
                                                           std::mutex & frictionConeMutex,
                                                           double m_frictionCoef,
                                                           double maxForce)
{
  unsigned int dim = 3;
  boost::shared_ptr<Polytope_Rn> newCone(new Polytope_Rn());
  // for now generate cone generates only the directions for the rays: we assume it is a polyhedral cone
  auto generators =
      generatePolyhedralConeGens(numberOfFrictionSides, contactSurface.rotation(), m_frictionCoef, maxForce);
  // here we manipulate polytope objects so need to add origin as a generator on the polyhedral cone
  generators.emplace_back(Eigen::Vector3d::Zero());
  for(const auto & g : generators)
  {
    boost::shared_ptr<Generator_Rn> gn(new Generator_Rn(dim));
    boost::numeric::ublas::vector<double> coords(3);
    coords.insert_element(0, g.x());
    coords.insert_element(1, g.y());
    coords.insert_element(2, g.z());
    gn->setCoordinates(coords);
    newCone->addGenerator(gn);
    // mc_rtc::log::info("Creating cone with vertex {}", g.transpose());
  }
  // update faces of the cone
  // The face computations are necessary for the minkowsky sum using normal fans
  DoubleDescriptionFromGenerators::Compute(newCone, 1000);
  // lock cone mutex, then reset cone pointer to newly computed cone
  std::lock_guard<std::mutex> lock(frictionConeMutex);
  frictionCone.reset();
  frictionCone = newCone;
  // mc_rtc::log::info("Created cone of dim {} with {} generators", forceCone->dimension(),
  //                   forceCone->numberOfGenerators());
}

// public
void DynamicPolytope::buildWrenchConeFromContact(int numberOfFrictionSides,
                                                 std::pair<std::pair<double, double>, sva::PTransformd> contactSurface,
                                                 boost::shared_ptr<Polytope_Rn> & forceCone,
                                                 boost::shared_ptr<Polytope_Rn> & momentPoly,
                                                 double m_frictionCoef,
                                                 double maxForce,
                                                 Eigen::Vector3d CoM)
{
  int dim = 3;
  boost::shared_ptr<Polytope_Rn> newForcePoly(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> newMomentPoly(new Polytope_Rn());
  // newCone
  // for now generate cone generates only the directions for the rays: we assume it is a polyhedral cone
  auto generators =
      generatePolyhedralConeGens(numberOfFrictionSides, contactSurface.second.rotation(), m_frictionCoef, maxForce);
  // contactPoint first is the pair of xHalfLength and yHalfLength of the rectangular contact
  std::vector<Eigen::Vector3d> points;
  points.emplace_back(contactSurface.first.first, contactSurface.first.second, 0);
  points.emplace_back(contactSurface.first.first, -contactSurface.first.second, 0);
  points.emplace_back(-contactSurface.first.first, -contactSurface.first.second, 0);
  points.emplace_back(-contactSurface.first.first, contactSurface.first.second, 0);

  // here we manipulate polytope objects so need to add origin as a generator on the polyhedral cone
  generators.emplace_back(Eigen::Vector3d::Zero());
  for(const auto & g : generators)
  {
    // this is the translational part: no variation of force depending on application point
    Eigen::Vector3d newForce = g;

    boost::shared_ptr<Generator_Rn> forceGN(new Generator_Rn(dim));
    boost::numeric::ublas::vector<double> forceCoords(3);
    forceCoords.insert_element(0, newForce.x());
    forceCoords.insert_element(1, newForce.y());
    forceCoords.insert_element(2, newForce.z());
    forceGN->setCoordinates(forceCoords);
    newForcePoly->addGenerator(forceGN);
    for(auto p : points)
    {
      // r is the extremity point of the surface : center + offset using half dimensions (-application point in world
      // to get desired frame)
      Eigen::Vector3d r = contactSurface.second.translation() + p - CoM;

      // here compute the generators using the "limits" of the surface contact (points r) and associate each generator
      // found for the contact this is the angular part of the 6d vector (resulting moment of the force generator at
      // application point, here CoM)
      Eigen::Vector3d newMoment = skewMatrix(r) * g;

      boost::shared_ptr<Generator_Rn> momentGN(new Generator_Rn(dim));
      boost::numeric::ublas::vector<double> momentCoords(3);

      momentCoords.insert_element(0, newMoment.x());
      momentCoords.insert_element(1, newMoment.y());
      momentCoords.insert_element(2, newMoment.z());
      momentGN->setCoordinates(momentCoords);
      newMomentPoly->addGenerator(momentGN);
    }
  }

  DoubleDescriptionFromGenerators::Compute(newForcePoly, 1000);
  DoubleDescriptionFromGenerators::Compute(newMomentPoly, 1000);

  forceCone->reset();
  forceCone = newForcePoly;
  momentPoly->reset();
  momentPoly = newMomentPoly;
}

// FIXME: restore implementation
void DynamicPolytope::computeFrictionConesFromContactSet(const mc_rbdyn::Robot & robot)
{
  // auto frictionCoeff = 0.7;
  // auto nbFrictionSides = 5;
  // auto maxForce = 250.;
  // auto CoM = robot.com();
  //
  // for(const auto & contactName : activeContacts_)
  // {
  //   // dummy value: should need difference between controlled surface and the target in the contact pair
  //   sva::PTransformd X_r1_r2 = sva::PTransformd::Identity();
  //   auto & lastContactPolytope = lastContactPolytopes_.at(contactName);
  //   if(!withMoments_)
  //   {
  //     buildFrictionConeFromContactWithHrep(nbFrictionSides, X_r1_r2, lastContactPolytope.frictionCone, frictionCoeff,
  //     3);
  //   }
  //   else
  //   {
  //     // find limits of contact area for moment limits
  //     double newContactHalfLength;
  //     double newContactHalfWidth;
  //     findHalfWidthLength(robot.surface(contactName), newContactHalfWidth, newContactHalfLength);
  //     std::pair<std::pair<double, double>, sva::PTransformd> newContact(
  //         std::pair<double, double>(newContactHalfLength, newContactHalfWidth), robot.surfacePose(contactName));
  //
  //     // TODO thread moments versions as well: add moment mutex + put mutexes as arguments
  //     buildWrenchConeFromContact(nbFrictionSides, newContact, lastContactPolytope.frictionCone,
  //         lastContactPolytope.frictionConeMoments, frictionCoeff, maxForce, CoM);
  //   }
  // }
}

//  TODO: restore implementation
void DynamicPolytope::computeForcePolyFromContactSet(const mc_rbdyn::Robot & robot)
{
  // double forceScalingFactor = 1;
  // for(const auto & contactName : activeContacts_)
  // {
  //   auto & lastContactPolytope = lastContactPolytopes_.at(contactName);
  //   buildActuationPolytopeFromContact(contactName, robot,
  //       lastContactPolytope.actuationPolytope, forceScalingFactor, 3);
  // }
}

void DynamicPolytope::computeFeasibleForcesFromContactSet(const mc_rbdyn::Robot & robot) {}

// TODO: restore implementation
void DynamicPolytope::computeMomentsRegion(Eigen::Vector3d /* comPosition */, const mc_rbdyn::Robot & robot)
{
  // // We scale the moment polytope according to the expression of the difference between eCMP and ZMP:
  // // eCMP = ZMP + 1/(m*(g+\ddot(c)_z)) * (tau_y, - tau_x, 0.)
  // // scale moment polytope by 1/(m*(g+\ddot(c)_z))
  // double scale = 1 / (robot.mass() * (9.81 + robot.comAcceleration().z()));
  // TopGeomTools::scalingFactor(CWCMoments_, scale);
  //
  // // TODO change coords from varignon (check) + check that corresponds to inside of eCMP region?
}

// Eigen::Vector3d projectPointInVRPRegion(Eigen::Vector3d testedPoint)
// {
//   constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(CWCForces_->getListOfGenerators());
//   for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
//     for (unsigned int j=0; j<iteGN.current()->numberOfFacets(); j++) {
//       boost::shared_ptr<HalfSpace_Rn> HS = iteGN.current()->getFacet(j);
//       boost::numeric::ublas::vector<double> projectedPoint;
//       double halfSpaceNorm = norm_2(HS->vect());
//       double disPoint2Hyp = HS->computeDistancePointHyperplane(iteGN.current()->vect(), projectedPoint,
//       halfSpaceNorm);
//       //std::cout << "d" << iteGN.currentIteratorNumber() << " = " << disPoint2Hyp << std::endl;
//       if (disPoint2Hyp > 0.25*TOL || disPoint2Hyp < -0.25*TOL)
//         isVeryClose = false;
//     }
//     if (isVeryClose == false) {
//       averagePoint /= iteGN.current()->numberOfFacets();
//       iteGN.current()->setCoordinates(averagePoint);
//     }
//   }
// }

// Eigen::Vector3d projectPointInPolytope(Eigen::Vector3d testedPoint, boost::shared_ptr<Polytope_Rn> & polytope)
// {
//   boost::numeric::ublas::vector<double> closestProjectedPoint;
//   constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(polytope->getListOfGenerators());
//   for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
//     for (unsigned int j=0; j<iteGN.current()->numberOfFacets(); j++) {
//       boost::shared_ptr<HalfSpace_Rn> HS = iteGN.current()->getFacet(j);
//       boost::numeric::ublas::vector<double> projectedPoint;
//       double halfSpaceNorm = norm_2(HS->vect());
//       double disPoint2Hyp = HS->computeDistancePointHyperplane(iteGN.current()->vect(), projectedPoint,
//       halfSpaceNorm);
//       //std::cout << "d" << iteGN.currentIteratorNumber() << " = " << disPoint2Hyp << std::endl;
//       if (disPoint2Hyp > 0.25*TOL || disPoint2Hyp < -0.25*TOL)
//         isVeryClose = false;
//     }
//     if (isVeryClose == false) {
//       averagePoint /= iteGN.current()->numberOfFacets();
//       iteGN.current()->setCoordinates(averagePoint);
//     }
//   }
//   return
// }

// TODO: restore implementation
bool DynamicPolytope::checkGravityCenterInPolytope(boost::shared_ptr<Polytope_Rn> & polytope)
{
  //   boost::numeric::ublas::vector<double> gravCenter(3);
  //   TopGeomTools::gravityCenter(polytope, gravCenter);
  //   Point_Rn testPoint(gravCenter(0), gravCenter(1), gravCenter(2));
  //
  //   int resultPolito = polytope->checkPoint(testPoint);
  //   bool retResultPolitopix = false;
  //   if(resultPolito == 1)
  //   {
  //     retResultPolitopix = true;
  //   }
  //
  //   Eigen::MatrixXd Normals;
  //   Eigen::VectorXd Offsets;
  //   updatePlanesMatrixConstraint(polytope, Normals, Offsets);
  //   Eigen::Vector3d testEigen(gravCenter(0), gravCenter(1), gravCenter(2));
  //   Eigen::VectorXd test = Normals * testEigen - Offsets;
  //   bool retResultEigen = true;
  //   for(int coeff = 0; coeff < test.size(); coeff++)
  //   {
  //     if(test(coeff) > 0.0)
  //     {
  //       retResultEigen = false;
  //     }
  //   }
  //
  //   return retResultPolitopix && retResultEigen;
}

// FIXME: should be done per-contact
// TODO: restore GUI computation
void DynamicPolytope::updateTrianglesContactsGUIPolitopix()
{
  // auto start_guiTriangles = mc_rtc::clock::now();
  // // Protecting set of names (threaded) then copying
  // contactSetMutex_.lock();
  // auto activeContactSet = activeContacts_;
  // auto contactsToBeRemoved = contactsToRemove_;
  // contactSetMutex_.unlock();
  //
  // int dim = Rn::getDimension() == 6 ? 6 : 3;
  //
  // for(const auto & contact : activeContactSet)
  // {
  //   auto contactPose = robot_.surfacePose(contact);
  //
  //   if(forcePolytopes_.at(contact)->dimension() == 3)
  //   {
  //     getContactMutex(frictionConeTrianglesMutexes_, contact).lock();
  //     getContactMutex(frictionConesMutexes_, contact).lock();
  //     update3DPolyTrianglesPolitopix(frictionCones_.at(contact), frictionConesTrianglesMap_.at(contact), guiScale_,
  //                                    contactPose);
  //     getContactMutex(frictionConesMutexes_, contact).unlock();
  //     getContactMutex(frictionConeTrianglesMutexes_, contact).unlock();
  //
  //     getContactMutex(forcePolyTrianglesMutexes_, contact).lock();
  //     getContactMutex(forcePolyMutexes_, contact).lock();
  //     update3DPolyTrianglesPolitopix(forcePolytopes_.at(contact), forcePolyTrianglesMap_.at(contact), guiScale_,
  //                                    contactPose);
  //     getContactMutex(forcePolyMutexes_, contact).unlock();
  //     getContactMutex(forcePolyTrianglesMutexes_, contact).unlock();
  //   }
  //   else
  //   {
  //     getContactMutex(momentTrianglesMutexes_, contact).lock();
  //     getContactMutex(forcePolyTrianglesMutexes_, contact).lock();
  //     update6DPolyTrianglesPolitopix(forcePolytopes_.at(contact), momentPolytopesTrianglesMap_.at(contact),
  //                                    forcePolyTrianglesMap_.at(contact), guiScale_, contactPose);
  //     getContactMutex(momentTrianglesMutexes_, contact).unlock();
  //     getContactMutex(forcePolyTrianglesMutexes_, contact).unlock();
  //   }
  // }
  // for(const auto & contact : contactsToBeRemoved)
  // {
  //   getContactMutex(frictionConeTrianglesMutexes_, contact).lock();
  //   frictionConesTrianglesMap_.at(contact).clear();
  //   getContactMutex(frictionConeTrianglesMutexes_, contact).unlock();
  //   if(dim == 6)
  //   {
  //     getContactMutex(momentTrianglesMutexes_, contact).lock();
  //     momentPolytopesTrianglesMap_.at(contact).clear();
  //     getContactMutex(momentTrianglesMutexes_, contact).unlock();
  //   }
  //   getContactMutex(forcePolyTrianglesMutexes_, contact).lock();
  //   forcePolyTrianglesMap_.at(contact).clear();
  //   getContactMutex(forcePolyTrianglesMutexes_, contact).unlock();
  // }
  // dt_compute_guiTrianglesContacts_ = mc_rtc::clock::now() - start_guiTriangles;
}

void DynamicPolytope::updateRBDynPolytopes(const Eigen::MatrixXd & Normals,
                                           const Eigen::VectorXd & Offsets,
                                           mc_rbdyn::Contact & contactRBDyn)
{
  mc_rbdyn::FeasiblePolytope polytope({Normals, Offsets});
  // Update the contact polytope of the desired robot
  if(robot_.robotIndex() == contactRBDyn.r1Index())
  {
    contactRBDyn.feasiblePolytopeR1(polytope);
  }
  else
  {
    contactRBDyn.feasiblePolytopeR2(polytope);
  }
}

// TODO: restore log
void DynamicPolytope::addToLogger(mc_rtc::Logger & logger)
{
  auto prefix = name_;
  logger_ = &logger;
  auto pp = "perf_" + prefix + "_MainThread_";
  logger.addLogEntry(pp + "totalLoop [ms]", this, [this]() { return dt_loop_total_.count(); });
  logger.addLogEntry(pp + "computeContactSet [ms]", this, [this]() { return dt_compute_contactSet_.count(); });
  logger.addLogEntry(pp + "dt_updateRBDynPolytopes [ms]", this, [this]() { return dt_update_rbdyn_.count(); });

  if(computeRegions_)
  {
    zmpRegionJob_.addToLogger(logger, name_ + "_ZMPRegion");
  }
  for(auto & [contactName, job] : feasiblePolytopesJobs_)
  {
    job.addToLogger(*logger_, name_ + "_" + contactName);
  }

  // }); logger.addLogEntry("perf_" + prefix + name_ + "_guiTrianglesContacts", this,
  //                    [this]() { return dt_guiTrianglesContacts().count(); });
  // logger.addLogEntry("perf_" + prefix + name_ + "_guiTrianglesRegions", this,
  //                    [this]() { return dt_guiTrianglesRegions().count(); });
}

void DynamicPolytope::addToGUI(mc_rtc::gui::StateBuilder * guiPtr)
{
  gui_ = guiPtr;

  auto & gui = *gui_;
  gui.removeCategory(guiCategory_);
  auto category = guiCategory_;

  if(computeRegions_)
  {
    auto zmpRegionCategory = category;
    zmpRegionCategory.push_back("ZMP Region");
    zmpRegionJob_.addToGUI(*gui_, guiCategory_);
  }

  auto coeffsCat = category;
  auto contactsCat = category;
  contactsCat.push_back("Contact Polytopes");
  auto CWCCat = category;
  CWCCat.push_back("Contact Wrench Cone");

  auto updateJobsGlobalParams = [this]()
  {
    for(auto & [_, job] : feasiblePolytopesJobs_)
    {
      job.HrepMode_ = HrepMode_;
      job.DDfrictionCones_ = DDfrictionCones_;
      job.combineWithFriction_ = combineWithFriction_;
    }
    VRPRegionMinkSumJob_.withMoments_ = withMoments_;
    VRPRegionMinkSumJob_.withVRPOffset_ = withVRPOffset_;
  };
  gui.addElement(this, category,
                 mc_rtc::gui::Checkbox(
                     "Compute explicit regions", [this]() { return computeRegions_; },
                     [this]()
                     {
                       computeRegions_ = !computeRegions_;
                       if(!computeRegions_)
                       {
                         VRPRegionMinkSumJob_.removeFromGUI();
                         zmpRegionJob_.removeFromGUI();
                       }
                     }),
                 mc_rtc::gui::Checkbox(
                     "Compute force poly from Hrep", [this]() { return HrepMode_; },
                     [this, updateJobsGlobalParams]()
                     {
                       HrepMode_ = !HrepMode_;
                       updateJobsGlobalParams();
                     }),
                 mc_rtc::gui::Checkbox(
                     "Compute moments", [this]() { return withMoments_; }, [this]() { withMoments_ = !withMoments_; }));

  // Debug options for the GUI
  gui.addElement(this, category,
                 mc_rtc::gui::Checkbox(
                     "Combine with frictions", [this]() { return combineWithFriction_; },
                     [this, updateJobsGlobalParams]()
                     {
                       combineWithFriction_ = !combineWithFriction_;
                       updateJobsGlobalParams();
                     }),
                 mc_rtc::gui::Checkbox(
                     "DD friction cones", [this]() { return DDfrictionCones_; },
                     [this, updateJobsGlobalParams]()
                     {
                       DDfrictionCones_ = !DDfrictionCones_;
                       updateJobsGlobalParams();
                     }),
                 mc_rtc::gui::Checkbox(
                     "With eCMP-VRP offset", [this]() { return withVRPOffset_; },
                     [this, updateJobsGlobalParams]()
                     {
                       withVRPOffset_ = !withVRPOffset_;
                       updateJobsGlobalParams();
                     }));
}

} // namespace mc_dynamic_polytopes
