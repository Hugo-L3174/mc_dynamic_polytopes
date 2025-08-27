#pragma once
#include <mc_rbdyn/Robot.h>
#include <SpaceVecAlg/SpaceVecAlg>
#include <future>
#include <optional>
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

  sva::PTransformd comPosition;
  std::vector<ZMPRegionContactInput> contactInputs;
};

struct ZMPRegionResult
{
  boost::shared_ptr<Polytope_Rn> zmpRegion;
};

struct ZMPRegionJob
{
  ZMPRegionInput input;

  void startAsync()
  {
    futureZMPRegion_ = std::async(std::bind(&ZMPRegionJob::computeZMPRegionJob, this));
    running_ = true;
  }

  bool checkResult()
  {
    if(futureZMPRegion_.wait_for(std::chrono::seconds(0)) == std::future_status::ready)
    {
      lastResult_ = futureZMPRegion_.get();
      running_ = false;
      return true;
    }
    return false;
  }

  /**
   * Get the last result if available
   */
  const std::optional<ZMPRegionResult> & lastResult() const noexcept
  {
    return lastResult_;
  }

  inline bool running() const noexcept
  {
    return running_;
  }

protected:
  bool running_ = false;
  std::future<ZMPRegionResult> futureZMPRegion_;
  std::optional<ZMPRegionResult> lastResult_;
  ZMPRegionResult computeZMPRegionJob();

protected: // actual computation
  boost::shared_ptr<Polytope_Rn> computeZMPRegion(Eigen::Vector3d comPosition);
};

} // namespace mc_dynamic_polytopes
