#define BOOST_TEST_MODULE DynamicPolytopeInstantiation
#include <mc_rbdyn/RobotLoader.h>
#include <mc_rtc/Configuration.h>
#include <boost/test/included/unit_test.hpp>
#include <mc_dynamic_polytopes/DynamicPolytope.h>

using namespace mc_dynamic_polytopes;

BOOST_AUTO_TEST_CASE(InstantiationTest)
{
  auto robots = mc_rbdyn::loadRobot(*mc_rbdyn::RobotLoader::get_robot_module("JVRC1"));
  auto & robot = robots->robot();
  mc_rtc::Configuration config;

  // Add required configuration keys if needed
  config.add("possibleContacts", std::vector<std::string>{"LeftFoot", "RightFoot"});
  config.add("withMoments", false);
  config.add("computeRegions", true);
  config.add("HrepMode", false);
  config.add("withFriction", true);

  BOOST_CHECK_NO_THROW(DynamicPolytope polytope("test", robot, config));
}
