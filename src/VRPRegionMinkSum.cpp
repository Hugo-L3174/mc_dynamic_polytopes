#include "mc_dynamic_polytopes/Time.h"
#include <mc_dynamic_polytopes/GUIComputations.h>
#include <mc_dynamic_polytopes/PolytopeFuncs.h>
#include <mc_dynamic_polytopes/VRPRegionMinkSum.h>

namespace mc_dynamic_polytopes
{

VRPRegionMinkSumJobResult VRPRegionMinkSumJob::computeJob()
{
  // XXX: anything to initialize in result_?
  // init CWC polytope
  result_.CWCForces.reset(new Polytope_Rn());
  result_.CWCMoments.reset(new Polytope_Rn());
  // init zmp region and intersection with ecmp region
  result_.zeroMomentRegion.reset(new Polytope_Rn());
  //
  // Step 3: All feasible regions computed, start minkowsky sum computation
  computeMinkowskySumPolitopix();

  // Step 4: wait for finished mink sum to convert to eCMP region (no thread needed)
  computeECMPRegion(input_.comPosition, input_.robotMass);

  // Step 5: wait for finished ZMP region to start zero moment intersection with eCMP region
  computeZeroMomentIntersection();

  // Step 6: translate eCMP region and zero-moment intersection to get VRP regions
  if(withVRPOffset_)
  {
    VRPtranslation(input_.comPosition.z());
  }

  // Update VRP planes internal variables to be fetched by controller
  auto start_updatePlanes = mc_rtc::clock::now();
  // VRPPlanesMutex_.lock();
  // CWCMutex_.lock();
  updatePlanesMatrixConstraint(result_.CWCForces, result_.DCMVRPPlanes.first, result_.DCMVRPPlanes.second);
  // CWCMutex_.unlock();
  // VRPPlanesMutex_.unlock();

  // Update zero moment region planes to be fetched by controller
  // zeroMomentPlanesMutex_.lock();
  // zeroMomentMutex_.lock();
  updatePlanesMatrixConstraint(result_.zeroMomentRegion, result_.zeroMomentPlanes.first,
                               result_.zeroMomentPlanes.second);
  // zeroMomentMutex_.unlock();
  // zeroMomentPlanesMutex_.unlock();
  vrpTimers_.dt_update_planes = mc_rtc::elapsed_ms_count(start_updatePlanes);

  // Update GUI display of regions
  updateTrianglesRegionsGUIPolitopix();
  return result_;
}

void VRPRegionMinkSumJob::computeMinkowskySumPolitopix()
{
  auto start_minkSum = mc_rtc::clock::now();
  boost::shared_ptr<Polytope_Rn> newForcePoly(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> newMomentPoly(new Polytope_Rn());

  // putting it in vector form for library function
  std::vector<boost::shared_ptr<Polytope_Rn>> polytopesForces;
  std::vector<boost::shared_ptr<Polytope_Rn>> polytopesMoments;

  // mc_rtc::log::info("Starting to add polytopes for mink sum");
  // for(const auto & active : activeContacts_)
  for(const auto & [active, contactPolytopes] : input_.contactsPolytopes)
  {
    boost::shared_ptr<Polytope_Rn> newContactPoly(new Polytope_Rn());
    // The friction + force intersection was overwritten in the force polytopes
    // Locking polytope mutex
    // Copying the contact polytopes to use for mink computation
    politopixAPI::copyPolytope(contactPolytopes.actuationPolytope, newContactPoly);
    // mc_rtc::log::info("Adding poly with {} gens and {} hs", newContactPoly->numberOfGenerators(),
    //                   newContactPoly->numberOfHalfSpaces());

    // After copying the contact frame polytope, rotate it to world frame before minkowski sum
    // Vectors are force in contact frame so X_contact_f
    // We want X_0_f = X_contact_f * X_0_contact
    auto X_0_contact = input_.contactsPose[active];
    X_0_contact.translation() = Eigen::Vector3d::Zero();

    // Rotating generators
    constIteratorOfListOfGeometricObjects<boost::shared_ptr<Generator_Rn>> iterGen(
        newContactPoly->getListOfGenerators());
    for(iterGen.begin(); iterGen.end() != true; iterGen.next())
    {
      Eigen::Vector3d vect(iterGen.current()->getCoordinate(0), iterGen.current()->getCoordinate(1),
                           iterGen.current()->getCoordinate(2));
      vect = (sva::PTransformd(vect) * X_0_contact).translation();
      iterGen.current()->setCoordinate(0, vect.x());
      iterGen.current()->setCoordinate(1, vect.y());
      iterGen.current()->setCoordinate(2, vect.z());
    }
    // Rotating halfspaces normals
    constIteratorOfListOfGeometricObjects<boost::shared_ptr<HalfSpace_Rn>> iterHS(
        newContactPoly->getListOfHalfSpaces());
    for(iterHS.begin(); iterHS.end() != true; iterHS.next())
    {
      Eigen::Vector3d normal(iterHS.current()->getCoefficient(0), iterHS.current()->getCoefficient(1),
                             iterHS.current()->getCoefficient(2));
      normal = (sva::PTransformd(normal) * X_0_contact).translation();
      iterHS.current()->setCoefficient(0, normal.x());
      iterHS.current()->setCoefficient(1, normal.y());
      iterHS.current()->setCoefficient(2, normal.z());
    }

    polytopesForces.emplace_back(newContactPoly);
    if(withMoments_)
    {
      polytopesMoments.emplace_back(contactPolytopes.frictionConeMoments);
    }
  }

  try
  {
    if(!polytopesForces.empty())
    {
      MinkowskiSum Mink(polytopesForces, newForcePoly);
    }
    result_.CWCForces = newForcePoly;
  }
  catch(const std::exception & e)
  {
    mc_rtc::log::error("[{}] Minkowski sum error: {}", name_, e.what());
  }

  // mc_rtc::log::info("CWCForces_ has {} generators and {} facets", CWCForces_->numberOfGenerators(),
  // CWCForces_->numberOfHalfSpaces());
  if(withMoments_)
  {
    MinkowskiSum Mink(polytopesMoments, newMomentPoly);
    // mc_rtc::log::info("CWCMoments_ has {} generators and {} facets", CWCMoments_->numberOfGenerators(),
    // CWCMoments_->numberOfHalfSpaces());
    result_.CWCMoments = newMomentPoly;
  }

  vrpTimers_.dt_compute_minkSum = mc_rtc::elapsed_ms_count(start_minkSum);
}

void VRPRegionMinkSumJob::computeECMPRegion(const Eigen::Vector3d & comPosition, double robotMass)
{
  /* We want to scale the force polytope according to the eCMP expression:
  eCMP = c - sumF/(m*(g/Dz)) --> eCMP = c - sumF*(Dz/mg)
  --> scale force polytope by - Deltaz/mg */

  // Several steps needed to mirror the polytope correctly (including the planes orientations):

  // First we scale with the POSITIVE scale because scaling only modifies the HS offsets, not their normals so the
  // polytope would be inverted but the planes inside out
  double scale = comPosition.z() / (robotMass * mc_rtc::constants::GRAVITY);
  TopGeomTools::scalingFactor(result_.CWCForces, scale);
  // Then we negate the polytope to mirror the generators and the normals correctly
  result_.CWCForces->negate();

  // Finally translate it from origin to the robot CoM to get eCMP region
  boost::numeric::ublas::vector<double> CoM(3);
  CoM[0] = comPosition.x();
  CoM[1] = comPosition.y();
  CoM[2] = comPosition.z();
  TopGeomTools::translate(result_.CWCForces, CoM);
}

void VRPRegionMinkSumJob::computeZeroMomentIntersection()
{
  auto start_zeroMomentIntersection = mc_rtc::clock::now();
  // making a deep copy of the force polytope to use as base for intersection with zmp region
  // (avoids shared_ptr problems)
  politopixAPI::copyPolytope(result_.CWCForces, result_.zeroMomentRegion);
  // politopixAPI::computeIntersection(CWCForces_, zmpRegion_, zeroMomentRegion_);
  politopixAPI::computeIntersectionWithoutCheck(result_.zeroMomentRegion, input_.zmpRegion);

  vrpTimers_.dt_zeroMoment_intersection = mc_rtc::elapsed_ms_count(start_zeroMomentIntersection);
}

void VRPRegionMinkSumJob::VRPtranslation(double deltaZ)
{
  boost::numeric::ublas::vector<double> deltaZVector(3);
  deltaZVector[0] = 0.;
  deltaZVector[1] = 0.;
  deltaZVector[2] = deltaZ;
  TopGeomTools::translate(result_.CWCForces, deltaZVector);
  TopGeomTools::translate(result_.zeroMomentRegion, deltaZVector);
}

void VRPRegionMinkSumJob::updateTrianglesRegionsGUIPolitopix()
{
  auto start_guiTriangles = mc_rtc::clock::now();
  auto noActiveContacts = input_.contactsPolytopes.empty(); // activeContacts_.empty();

  if(!noActiveContacts)
  {
    // gui scale for CWC should be 1, it is position space and not force space (because eCMP)
    // update6DPolyTrianglesPolitopix(CWC_, CWCMomentTriangles_, CWCForceTriangles_, guiScale_);
    update3DPolyTrianglesPolitopix(result_.CWCForces, result_.CWCForceTriangles, 1);

    // scale 1 here: already position space
    update3DPolyTrianglesPolitopix(input_.zmpRegion, result_.ZMPTriangles, 1);
    update3DPolyTrianglesPolitopix(result_.zeroMomentRegion, result_.zeroMomentTriangles, 1);

    if(withMoments_)
    {
      update3DPolyTrianglesPolitopix(result_.CWCMoments, result_.CWCMomentTriangles, guiScale_);
    }
  }
  else
  {
    result_.CWCForceTriangles.clear();

    result_.ZMPTriangles.clear();

    result_.zeroMomentTriangles.clear();

    if(withMoments_)
    {
      result_.CWCMomentTriangles.clear();
    }
  }
  vrpTimers_.dt_compute_guiTrianglesRegions_ = mc_rtc::elapsed_ms_count(start_guiTriangles);
}

} // namespace mc_dynamic_polytopes
