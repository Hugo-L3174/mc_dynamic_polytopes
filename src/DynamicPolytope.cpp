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
  possibleContacts_ = dynamicPolyConfig("possibleContacts", std::set<std::string>{"LeftFoot", "RightFoot"});
  withMoments_ = dynamicPolyConfig("withMoments", false);
  computeRegions_ = dynamicPolyConfig("computeRegions", true);
  HrepMode_ = dynamicPolyConfig("HrepMode", false);
  combineWithFriction_ = dynamicPolyConfig("withFriction", true);

  // Init dimension
  unsigned int dim = withMoments_ ? 6 : 3;
  Rn::setDimension(dim);
  Rn::setTolerance(1.e-07);

  for(const auto & contact : possibleContacts_)
  {
    // init triangles maps
    std::vector<std::array<Eigen::Vector3d, 3UL>> newTrianglesArray;
    frictionConesTrianglesMap_.emplace(contact, newTrianglesArray);
    forcePolyTrianglesMap_.emplace(contact, newTrianglesArray);
    momentPolytopesTrianglesMap_.emplace(contact, newTrianglesArray);

    refContactTransforms_.emplace(contact, sva::PTransformd::Identity());
  }
  // init CWC polytope
  CWCForces_.reset(new Polytope_Rn());
  CWCMoments_.reset(new Polytope_Rn());

  // init zmp region and intersection with ecmp region
  zmpRegion_.reset(new Polytope_Rn());
  zeroMomentRegion_.reset(new Polytope_Rn());
}

DynamicPolytope::~DynamicPolytope()
{
  stopThread();
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
    job.load(config_); // XXX should load per-contact config instead of global one
    job.contactTimers.dt_constructor = mc_rtc::elapsed_ms_count(start_constructor);
    job.HrepMode_ = HrepMode_; // XXX

    if(!job.running())
    {
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
                                              // XXX: this uses the default value anyway
                                              // input.numberOfFrictionSides = numbersOfFrictionSides_.at(contactName);
                                              // input.forceScalingFactor = forceScalingFactors_.at(contactName);
                                              input.frictionCoefficient = contactsRBDyn_.at(contactName).friction();
                                            })
                                            .count();

      // Job starts as an async task, use job.checkResult() later to know whether it is finished and retrive its value
      job.startAsync();
      if(logger_)
      {
        job.addToLogger(*logger_, name_);
      }
      if(gui_)
      {
        job.addToGUI(*gui_, guiCategory_);
      }
    }
  }

  // Step 2: wait for finished individual feasible regions
  // Wait for finished threads and join them, then update their matrix constraints
  auto start_debug = mc_rtc::clock::now();
  for(const auto & contactName : activeContacts_)
  {
    auto & job = feasiblePolytopesJobs_[contactName];
    if(job.checkResult())
    {
      const auto & result = *job.lastResult();
      updateRBDynPolytopes(result.forcePolyPlanes.first, result.forcePolyPlanes.second, contactsRBDyn_.at(contactName));
    }
  }
  dt_debug_ = mc_rtc::elapsed_ms(start_debug);

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
    }
  }

  dt_compute_contactSet_ = mc_rtc::elapsed_ms(start_loop);

  // Steps 3-6: launch the rest everytime the previous full region was computed
  // if(computeRegions_) // TODO: restore implementation
  // {
  //   if(VRPRegionComputed_)
  //   {
  //     if(minkSumThread_.joinable())
  //     {
  //       minkSumThread_.join();
  //     }
  //
  //     VRPRegionComputed_ = false;
  //     minkSumThread_ = std::thread(&DynamicPolytope::computeVRPRegionWithMinkSum, this);
  //   }
  // }

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
  for(auto g : generators)
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
void DynamicPolytope::computeMinkowskySumPolitopix()
{
  // auto start_minkSum = mc_rtc::clock::now();
  // boost::shared_ptr<Polytope_Rn> newForcePoly(new Polytope_Rn());
  // boost::shared_ptr<Polytope_Rn> newMomentPoly(new Polytope_Rn());
  //
  // // putting it in vector form for library function
  // std::vector<boost::shared_ptr<Polytope_Rn>> polytopesForces;
  // std::vector<boost::shared_ptr<Polytope_Rn>> polytopesMoments;
  //
  // // protecting set of contact names
  // contactSetMutex_.lock();
  // // mc_rtc::log::info("Starting to add polytopes for mink sum");
  // for(const auto & active : activeContacts_)
  // {
  //   boost::shared_ptr<Polytope_Rn> newContactPoly(new Polytope_Rn());
  //   // The friction + force intersection was overwritten in the force polytopes
  //   // Locking polytope mutex
  //   // Copying the contact polytopes to use for mink computation
  //   getContactMutex(forcePolyMutexes_, active).lock();
  //   politopixAPI::copyPolytope(forcePolytopes_.at(active), newContactPoly);
  //   // mc_rtc::log::info("Adding poly with {} gens and {} hs", newContactPoly->numberOfGenerators(),
  //   //                   newContactPoly->numberOfHalfSpaces());
  //   getContactMutex(forcePolyMutexes_, active).unlock();
  //
  //   // After copying the contact frame polytope, rotate it to world frame before minkowski sum
  //   // Vectors are force in contact frame so X_contact_f
  //   // We want X_0_f = X_contact_f * X_0_contact
  //   auto X_0_contact = robot_.surfacePose(active);
  //   X_0_contact.translation() = Eigen::Vector3d::Zero();
  //
  //   // Rotating generators
  //   constIteratorOfListOfGeometricObjects<boost::shared_ptr<Generator_Rn>> iterGen(
  //       newContactPoly->getListOfGenerators());
  //   for(iterGen.begin(); iterGen.end() != true; iterGen.next())
  //   {
  //     Eigen::Vector3d vect(iterGen.current()->getCoordinate(0), iterGen.current()->getCoordinate(1),
  //                          iterGen.current()->getCoordinate(2));
  //     vect = (sva::PTransformd(vect) * X_0_contact).translation();
  //     iterGen.current()->setCoordinate(0, vect.x());
  //     iterGen.current()->setCoordinate(1, vect.y());
  //     iterGen.current()->setCoordinate(2, vect.z());
  //   }
  //   // Rotating halfspaces normals
  //   constIteratorOfListOfGeometricObjects<boost::shared_ptr<HalfSpace_Rn>> iterHS(
  //       newContactPoly->getListOfHalfSpaces());
  //   for(iterHS.begin(); iterHS.end() != true; iterHS.next())
  //   {
  //     Eigen::Vector3d normal(iterHS.current()->getCoefficient(0), iterHS.current()->getCoefficient(1),
  //                            iterHS.current()->getCoefficient(2));
  //     normal = (sva::PTransformd(normal) * X_0_contact).translation();
  //     iterHS.current()->setCoefficient(0, normal.x());
  //     iterHS.current()->setCoefficient(1, normal.y());
  //     iterHS.current()->setCoefficient(2, normal.z());
  //   }
  //
  //   polytopesForces.emplace_back(newContactPoly);
  //   if(withMoments_)
  //   {
  //     polytopesMoments.emplace_back(frictionConesMoments_.at(active));
  //   }
  // }
  // contactSetMutex_.unlock();
  //
  // try
  // {
  //   if(!polytopesForces.empty())
  //   {
  //     MinkowskiSum Mink(polytopesForces, newForcePoly);
  //   }
  //   std::lock_guard<std::mutex> lock(CWCMutex_);
  //   CWCForces_.reset();
  //   CWCForces_ = newForcePoly;
  // }
  // catch(const std::exception & e)
  // {
  //   mc_rtc::log::error("[{}] Minkowski sum error: {}", name_, e.what());
  // }
  //
  // // mc_rtc::log::info("CWCForces_ has {} generators and {} facets", CWCForces_->numberOfGenerators(),
  // // CWCForces_->numberOfHalfSpaces());
  // if(withMoments_)
  // {
  //   MinkowskiSum Mink(polytopesMoments, newMomentPoly);
  //   // mc_rtc::log::info("CWCMoments_ has {} generators and {} facets", CWCMoments_->numberOfGenerators(),
  //   // CWCMoments_->numberOfHalfSpaces());
  //   CWCMoments_.reset();
  //   CWCMoments_ = newMomentPoly;
  // }
  //
  // dt_compute_minkSum_ = mc_rtc::clock::now() - start_minkSum;
}

void DynamicPolytope::computeECMPRegion(Eigen::Vector3d comPosition, const mc_rbdyn::Robot & robot)
{
  /* We want to scale the force polytope according to the eCMP expression:
  eCMP = c - sumF/(m*(g/Dz)) --> eCMP = c - sumF*(Dz/mg)
  --> scale force polytope by - Deltaz/mg */

  // Several steps needed to mirror the polytope correctly (including the planes orientations):

  // First we scale with the POSITIVE scale because scaling only modifies the HS offsets, not their normals so the
  // polytope would be inverted but the planes inside out
  double scale = comPosition.z() / (robot.mass() * 9.81);
  TopGeomTools::scalingFactor(CWCForces_, scale);
  // Then we negate the polytope to mirror the generators and the normals correctly
  CWCForces_->negate();

  // Finally translate it from origin to the robot CoM to get eCMP region
  boost::numeric::ublas::vector<double> CoM(3);
  CoM[0] = comPosition.x();
  CoM[1] = comPosition.y();
  CoM[2] = comPosition.z();
  TopGeomTools::translate(CWCForces_, CoM);
}

void DynamicPolytope::VRPtranslation(double deltaZ)
{
  boost::numeric::ublas::vector<double> deltaZVector(3);
  deltaZVector[0] = 0.;
  deltaZVector[1] = 0.;
  deltaZVector[2] = deltaZ;
  std::lock_guard<std::mutex> lockCWC(CWCMutex_);
  TopGeomTools::translate(CWCForces_, deltaZVector);
  std::lock_guard<std::mutex> lockZeroMoment(zeroMomentMutex_);
  TopGeomTools::translate(zeroMomentRegion_, deltaZVector);
}

void DynamicPolytope::computeZeroMomentIntersection()
{
  auto start_zeroMomentIntersection = mc_rtc::clock::now();
  zeroMomentMutex_.lock();
  zeroMomentRegion_->reset();
  // making a deep copy of the force polytope to use as base for intersection with zmp region
  // (avoids shared_ptr problems)
  CWCMutex_.lock();
  politopixAPI::copyPolytope(CWCForces_, zeroMomentRegion_);
  CWCMutex_.unlock();

  ZMPMutex_.lock();
  // politopixAPI::computeIntersection(CWCForces_, zmpRegion_, zeroMomentRegion_);
  politopixAPI::computeIntersectionWithoutCheck(zeroMomentRegion_, zmpRegion_);
  ZMPMutex_.unlock();
  zeroMomentMutex_.unlock();

  dt_zeroMoment_intersection_ = mc_rtc::clock::now() - start_zeroMomentIntersection;
}

// TODO: restore implementation
void DynamicPolytope::computeVRPRegionWithMinkSum()
{
  // // Step 3: All feasible regions computed, start minkowsky sum computation
  // computeMinkowskySumPolitopix();
  //
  // // Step 4: wait for finished mink sum to convert to eCMP region (no thread needed)
  // computeECMPRegion(robot_.com(), robot_);
  //
  // // Step 5: wait for finished ZMP region to start zero moment intersection with eCMP region
  // if(zmpThread_.joinable())
  // {
  //   zmpThread_.join();
  // }
  // computeZeroMomentIntersection();
  //
  // // Step 6: translate eCMP region and zero-moment intersection to get VRP regions
  // if(withVRPOffset_)
  // {
  //   VRPtranslation(robot_.com().z());
  // }
  //
  // // Update VRP planes internal variables to be fetched by controller
  // auto start_updatePlanes = mc_rtc::clock::now();
  // VRPPlanesMutex_.lock();
  // CWCMutex_.lock();
  // updatePlanesMatrixConstraint(CWCForces_, DCMVRPPlanes_.first, DCMVRPPlanes_.second);
  // CWCMutex_.unlock();
  // VRPPlanesMutex_.unlock();
  //
  // // Update zero moment region planes to be fetched by controller
  // zeroMomentPlanesMutex_.lock();
  // zeroMomentMutex_.lock();
  // updatePlanesMatrixConstraint(zeroMomentRegion_, zeroMomentPlanes_.first, zeroMomentPlanes_.second);
  // zeroMomentMutex_.unlock();
  // zeroMomentPlanesMutex_.unlock();
  // dt_update_planes_ = mc_rtc::clock::now() - start_updatePlanes;
  //
  // // Update GUI display of regions
  // updateTrianglesRegionsGUIPolitopix();
  // VRPRegionComputed_ = true;
}

void DynamicPolytope::computeMomentsRegion(Eigen::Vector3d /* comPosition */, const mc_rbdyn::Robot & robot)
{
  // We scale the moment polytope according to the expression of the difference between eCMP and ZMP:
  // eCMP = ZMP + 1/(m*(g+\ddot(c)_z)) * (tau_y, - tau_x, 0.)
  // scale moment polytope by 1/(m*(g+\ddot(c)_z))
  double scale = 1 / (robot.mass() * (9.81 + robot.comAcceleration().z()));
  TopGeomTools::scalingFactor(CWCMoments_, scale);

  // TODO change coords from varignon (check) + check that corresponds to inside of eCMP region?
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

void DynamicPolytope::updateTrianglesRegionsGUIPolitopix()
{
  auto start_guiTriangles = mc_rtc::clock::now();
  contactSetMutex_.lock();
  auto noActiveContacts = activeContacts_.empty();
  contactSetMutex_.unlock();

  if(!noActiveContacts)
  {
    // gui scale for CWC should be 1, it is position space and not force space (because eCMP)
    // update6DPolyTrianglesPolitopix(CWC_, CWCMomentTriangles_, CWCForceTriangles_, guiScale_);
    CWCForceTrianglesMutex_.lock();
    CWCMutex_.lock();
    update3DPolyTrianglesPolitopix(CWCForces_, CWCForceTriangles_, 1);
    CWCMutex_.unlock();
    CWCForceTrianglesMutex_.unlock();

    // scale 1 here: already position space
    ZMPTrianglesMutex_.lock();
    ZMPMutex_.lock();
    update3DPolyTrianglesPolitopix(zmpRegion_, ZMPTriangles_, 1);
    ZMPMutex_.unlock();
    ZMPTrianglesMutex_.unlock();

    zeroMomentTrianglesMutex_.lock();
    zeroMomentMutex_.lock();
    update3DPolyTrianglesPolitopix(zeroMomentRegion_, zeroMomentTriangles_, 1);
    zeroMomentMutex_.unlock();
    zeroMomentTrianglesMutex_.unlock();

    if(withMoments_)
    {
      CWCMomentTrianglesMutex_.lock();
      update3DPolyTrianglesPolitopix(CWCMoments_, CWCMomentTriangles_, guiScale_);
      CWCMomentTrianglesMutex_.unlock();
    }
  }
  else
  {
    CWCForceTrianglesMutex_.lock();
    CWCForceTriangles_.clear();
    CWCForceTrianglesMutex_.unlock();

    ZMPTrianglesMutex_.lock();
    ZMPTriangles_.clear();
    ZMPTrianglesMutex_.unlock();

    zeroMomentTrianglesMutex_.lock();
    zeroMomentTriangles_.clear();
    zeroMomentTrianglesMutex_.unlock();

    if(withMoments_)
    {
      CWCMomentTrianglesMutex_.lock();
      CWCMomentTriangles_.clear();
      CWCMomentTrianglesMutex_.unlock();
    }
  }
  dt_compute_guiTrianglesRegions_ = mc_rtc::clock::now() - start_guiTriangles;
}

void DynamicPolytope::updateRBDynPolytopes(const Eigen::MatrixXd & Normals,
                                           const Eigen::VectorXd & Offsets,
                                           mc_rbdyn::Contact & contactRBDyn)
{
  mc_rbdyn::FeasiblePolytope polytope({Normals, Offsets});
  // Update the contact polytope of the desired robot
  // FIXME : This contact can be dead at this point, find a way to catch the segfault
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
  auto pp = "perf_" + prefix + "_";
  logger.addLogEntry(pp + "totalLoop [ms]", this, [this]() { return dt_loop_total_.count(); });
  logger.addLogEntry(pp + "computeContactSet [ms]", this, [this]() { return dt_compute_contactSet_.count(); });
  logger.addLogEntry(pp + "dt_debug [ms]", this, [this]() { return dt_debug_.count(); });

  for(auto & [_, job] : feasiblePolytopesJobs_)
  {
    job.addToLogger(*logger_, name_);
  }

  // logger.addLogEntry("perf_" + prefix + name_ + "_minkSum", this, [this]() { return dt_minkSum().count(); });
  // logger.addLogEntry("perf_" + prefix + name_ + "_updatePlanes", this, [this]() { return dt_updatePlanes().count();
  // }); logger.addLogEntry("perf_" + prefix + name_ + "_guiTrianglesContacts", this,
  //                    [this]() { return dt_guiTrianglesContacts().count(); });
  // logger.addLogEntry("perf_" + prefix + name_ + "_guiTrianglesRegions", this,
  //                    [this]() { return dt_guiTrianglesRegions().count(); });
  // logger.addLogEntry("perf_" + prefix + name_ + "_zeroMomentIntersection", this,
  //                    [this]() { return dt_zeroMomentInter().count(); });
  // for(const auto & contact : possibleContacts_)
  // {
  //   logger.addLogEntry("perf_" + prefix + name_ + "_" + contact + "_frictionCone", this,
  //                      [this, contact]() { return getContact_dt_frictionCone(contact).count(); });
  //   logger.addLogEntry("perf_" + prefix + name_ + "_" + contact + "_forcePolytope", this,
  //                      [this, contact]() { return getContact_dt_forcePolytope(contact).count(); });
  //   logger.addLogEntry("perf_" + prefix + name_ + "_" + contact + "_intersection", this,
  //                      [this, contact]() { return getContact_dt_intersection(contact).count(); });
  //   logger.addLogEntry("perf_" + prefix + name_ + "_" + contact + "_Total", this,
  //                      [this, contact]() { return getContact_dt_Total(contact).count(); });
  // }
}

void DynamicPolytope::addToGUI(mc_rtc::gui::StateBuilder * guiPtr)
{
  gui_ = guiPtr;

  auto & gui = *gui_;
  gui.removeCategory(guiCategory_);
  auto category = guiCategory_;

  // XXX: simplify this
  config_("gui")("polyhedronForce")("triangle_color", polyForceConfig_.triangle_color);
  config_("gui")("polyhedronForce")("show_triangle", polyForceConfig_.show_triangle);
  config_("gui")("polyhedronForce")("use_triangle_color", polyForceConfig_.use_triangle_color);
  config_("gui")("polyhedronForce")("edges", polyForceConfig_.edge_config);
  config_("gui")("polyhedronForce")("show_edges", polyForceConfig_.show_edges);
  config_("gui")("polyhedronForce")("fixed_edge_color", polyForceConfig_.fixed_edge_color);
  config_("gui")("polyhedronForce")("vertices")("color", polyForceConfig_.vertices_config.color);
  config_("gui")("polyhedronForce")("vertices")("scale", polyForceConfig_.vertices_config.scale);
  config_("gui")("polyhedronForce")("show_vertices", polyForceConfig_.show_vertices);
  config_("gui")("polyhedronForce")("fixed_vertices_color", polyForceConfig_.fixed_vertices_color);

  config_("gui")("polyhedronMoment")("triangle_color", polyMomentConfig_.triangle_color);
  config_("gui")("polyhedronMoment")("show_triangle", polyMomentConfig_.show_triangle);
  config_("gui")("polyhedronMoment")("use_triangle_color", polyMomentConfig_.use_triangle_color);
  config_("gui")("polyhedronMoment")("edges", polyMomentConfig_.edge_config);
  config_("gui")("polyhedronMoment")("show_edges", polyMomentConfig_.show_edges);
  config_("gui")("polyhedronMoment")("fixed_edge_color", polyMomentConfig_.fixed_edge_color);
  config_("gui")("polyhedronMoment")("vertices")("color", polyMomentConfig_.vertices_config.color);
  config_("gui")("polyhedronMoment")("vertices")("scale", polyMomentConfig_.vertices_config.scale);
  config_("gui")("polyhedronMoment")("show_vertices", polyMomentConfig_.show_vertices);
  config_("gui")("polyhedronMoment")("fixed_vertices_color", polyMomentConfig_.fixed_vertices_color);

  config_("gui")("polyhedronZMP")("triangle_color", polyZMPConfig_.triangle_color);
  config_("gui")("polyhedronZMP")("show_triangle", polyZMPConfig_.show_triangle);
  config_("gui")("polyhedronZMP")("use_triangle_color", polyZMPConfig_.use_triangle_color);
  config_("gui")("polyhedronZMP")("edges", polyZMPConfig_.edge_config);
  config_("gui")("polyhedronZMP")("show_edges", polyZMPConfig_.show_edges);
  config_("gui")("polyhedronZMP")("fixed_edge_color", polyZMPConfig_.fixed_edge_color);
  config_("gui")("polyhedronZMP")("vertices")("color", polyZMPConfig_.vertices_config.color);
  config_("gui")("polyhedronZMP")("vertices")("scale", polyZMPConfig_.vertices_config.scale);
  config_("gui")("polyhedronZMP")("show_vertices", polyZMPConfig_.show_vertices);
  config_("gui")("polyhedronZMP")("fixed_vertices_color", polyZMPConfig_.fixed_vertices_color);

  config_("gui")("polyhedronZeroMomentArea")("triangle_color", polyZeroMomentAreaConfig_.triangle_color);
  config_("gui")("polyhedronZeroMomentArea")("show_triangle", polyZeroMomentAreaConfig_.show_triangle);
  config_("gui")("polyhedronZeroMomentArea")("use_triangle_color", polyZeroMomentAreaConfig_.use_triangle_color);
  config_("gui")("polyhedronZeroMomentArea")("edges", polyZeroMomentAreaConfig_.edge_config);
  config_("gui")("polyhedronZeroMomentArea")("show_edges", polyZeroMomentAreaConfig_.show_edges);
  config_("gui")("polyhedronZeroMomentArea")("fixed_edge_color", polyZeroMomentAreaConfig_.fixed_edge_color);
  config_("gui")("polyhedronZeroMomentArea")("vertices")("color", polyZeroMomentAreaConfig_.vertices_config.color);
  config_("gui")("polyhedronZeroMomentArea")("vertices")("scale", polyZeroMomentAreaConfig_.vertices_config.scale);
  config_("gui")("polyhedronZeroMomentArea")("show_vertices", polyZeroMomentAreaConfig_.show_vertices);
  config_("gui")("polyhedronZeroMomentArea")("fixed_vertices_color", polyZeroMomentAreaConfig_.fixed_vertices_color);

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
      job.guiScale_ = guiScale_;
      job.combineWithFriction_ = combineWithFriction_;
    }
  };
  gui.addElement(this, category,
                 mc_rtc::gui::Checkbox(
                     "Compute explicit regions", [this]() { return computeRegions_; },
                     [this]() { computeRegions_ = !computeRegions_; }),
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
                     }));
  // mc_rtc::gui::Checkbox(
  //     "With eCMP-VRP offset", [this]() { return withVRPOffset_; }, [this]() { withVRPOffset_ = !withVRPOffset_;
  //     }));
  //
  // gui.addElement(
  //     this, CWCCat,
  //     mc_rtc::gui::Polyhedron("CWC forces", polyForceConfig_, [this]() { return getCWCForceTriangles(); }),
  //     /*mc_rtc::gui::Polyhedron("CWC moments", polyMomentConfig_, [this]() { return getCWCMomentTriangles(); }),*/
  //     mc_rtc::gui::Polyhedron("ZMP area", polyZMPConfig_, [this]() { return getZMPTriangles(); }),
  //     mc_rtc::gui::Polyhedron("Zero moment region", polyZeroMomentAreaConfig_,
  //                             [this]() { return getZeroMomentTriangles(); }));
}

} // namespace mc_dynamic_polytopes
