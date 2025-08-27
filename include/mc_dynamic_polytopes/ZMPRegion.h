#pragma once
#include <mc_rbdyn/Robot.h>
#include <SpaceVecAlg/SpaceVecAlg>
#include <mc_dynamic_polytopes/AsyncJob.h>
#include <politopix/politopixAPI.h>

namespace mc_dynamic_polytopes
{

struct ZMPRegionContactInput
{
  std::string contactName;
  std::vector<sva::PTransformd> points;
  sva::PTransformd X_0_s;
  sva::PTransformd X_b_s;
};

struct ZMPRegionInput
{
  ZMPRegionInput() = default;

  void initialize(const mc_rbdyn::Robot & robot, const std::set<std::string> & contacts)
  {
    comPosition = robot.com();
    for(const auto & contactName : contacts)
    {
      if(!robot.hasSurface(contactName))
      {
        mc_rtc::log::error_and_throw("ZMPRegionInput: robot has no surface named {}", contactName);
      }
      const auto & surface = robot.surface(contactName);
      ZMPRegionContactInput contactInput;
      contactInput.contactName = contactName;
      contactInput.points = surface.points();
      contactInput.X_0_s = robot.surfacePose(contactName);
      contactInput.X_b_s = surface.X_b_s();
      contactInputs.push_back(contactInput);
    }
  }

  Eigen::Vector3d comPosition;
  std::vector<ZMPRegionContactInput> contactInputs;
};

struct ZMPRegionResult
{
  boost::shared_ptr<Polytope_Rn> zmpRegion;
};

struct ZMPRegionJob : MakeAsyncJob<ZMPRegionJob, ZMPRegionInput, ZMPRegionResult>
{
  // Implements the required CRTP method
  ZMPRegionResult computeJob()
  {
    ZMPRegionResult result;
    result.zmpRegion = computeZMPRegion(input.comPosition);
    return result;
  }

  // Optionally implement logger/GUI integration
  void addToLoggerImpl() {}
  void addToGUIImpl() {}

protected:
  boost::shared_ptr<Polytope_Rn> computeZMPRegion(Eigen::Vector3d comPosition);
};

} // namespace mc_dynamic_polytopes
