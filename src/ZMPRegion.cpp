#include <mc_dynamic_polytopes/PolytopeFuncs.h>
#include <mc_dynamic_polytopes/ZMPRegion.h>

namespace mc_dynamic_polytopes
{

ZMPRegionResult ZMPRegionJob::computeZMPRegionJob()
{
  ZMPRegionResult result;
  result.zmpRegion = computeZMPRegion(input.comPosition.translation());
  return result;
}

boost::shared_ptr<Polytope_Rn> ZMPRegionJob::computeZMPRegion(Eigen::Vector3d comPosition)
{
  // XXX dummy zone for now: convex area formed by the polygon envelope of feet + com position
  int dim = 3;
  boost::shared_ptr<Polytope_Rn> newZMPPoly(new Polytope_Rn());

  std::vector<double> qhullVect;
  const std::string contactTest = "LeftFoot";
  for(const auto & contactInput : input.contactInputs)
  {
    // Important ! points for surfaces give the points coordinates from the parent link, not from the
    // surface origin, so X_b_p and not X_s_p
    // This means to get world frame surface points X_0_p = X_s_p * X_0_s we need X_s_p
    // X_s_p = X_b_p * X_s_b = X_b_p * X_b_s.inv()
    const auto & cPoints = contactInput.points;
    for(auto cPoint : cPoints)
    {
      cPoint = cPoint * contactInput.X_b_s.inv() * contactInput.X_0_s;
      // coords.insert_element(0, cPoint.translation().x());
      // coords.insert_element(1, cPoint.translation().y());
      // coords.insert_element(2, cPoint.translation().z());
      // Testing with triple distance points (they are generators for a polyhedral cone, so they should not be bounded)
      qhullVect.emplace_back(3 * cPoint.translation().x() - 2 * comPosition.x());
      qhullVect.emplace_back(3 * cPoint.translation().y() - 2 * comPosition.y());
      qhullVect.emplace_back(3 * cPoint.translation().z() - 2 * comPosition.z());
    }
  }

  // Use this to test only one contact
  /*
  auto cPoints = robot_.surface(contactTest).points();
  for(auto cPoint : cPoints)
  {
    cPoint = cPoint * robot_.surface(contactTest).X_b_s().inv() * robot_.surface(contactTest).X_0_s(robot_);
    boost::shared_ptr<Generator_Rn> gn(new Generator_Rn(dim));
    boost::numeric::ublas::vector<double> coords(3);
    coords.insert_element(0, cPoint.translation().x());
    coords.insert_element(1, cPoint.translation().y());
    coords.insert_element(2, cPoint.translation().z());
    gn->setCoordinates(coords);
    newZMPPoly->addGenerator(gn);
  }
  */

  qhullVect.emplace_back(comPosition.x());
  qhullVect.emplace_back(comPosition.y());
  qhullVect.emplace_back(comPosition.z());

  computeQhullHrep(qhullVect, newZMPPoly, dim);

  politopixAPI::computeDoubleDescriptionWithoutCheck(newZMPPoly, 3000);

  checkAllHSInternal("ZMP", newZMPPoly);

  return newZMPPoly;
}

} // namespace mc_dynamic_polytopes
