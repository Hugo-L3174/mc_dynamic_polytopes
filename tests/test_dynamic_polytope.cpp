#define BOOST_TEST_MODULE DynamicPolytopeInstantiation
#include <mc_rbdyn/RobotLoader.h>
#include <mc_rtc/Configuration.h>
#include <boost/test/included/unit_test.hpp>
#include <DynamicPolytope.h>

BOOST_AUTO_TEST_CASE(InstantiationTest)
{
  // Minimal robot and configuration setup
  // Replace with actual robot loading if available
  // mc_rbdyn::Robot robot("test_robot", mc_rbdyn::RobotModulePtr{});
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
