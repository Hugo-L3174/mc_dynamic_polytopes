/* @todo
- Try out fake robots instanciations with chosen number of dofs like in mc_rtc unit tests
- Check if total number of dofs in the robot has an impact (size of jacobian and inertia matrices)
*/

#include <mc_rbdyn/RobotLoader.h>
#include <mc_rtc/Configuration.h>
#include <benchmark/benchmark.h>
#include <mc_dynamic_polytopes/DynamicPolytope.h>

static void BM_PolytopeJVRC1_3D_6DOF(benchmark::State & state)
{
  // Setup part (not timed)
  auto robots = mc_rbdyn::loadRobot(*mc_rbdyn::RobotLoader::get_robot_module("JVRC1"));
  auto & robot = robots->robot();
  auto envs = mc_rbdyn::loadRobot(*mc_rbdyn::RobotLoader::get_robot_module("env/ground"));
  auto & ground = envs->robot();

  mc_rbdyn::RobotsPtr robotsList(mc_rbdyn::Robots::make());
  robotsList->robotCopy(robot, robot.name());
  robotsList->robotCopy(ground, ground.name());

  mc_rtc::Configuration config;
  // Add required configuration keys if needed
  config.add("possibleContacts", std::vector<std::string>{"LeftFoot", "RightFoot"});
  config.add("withMoments", false);
  config.add("computeRegions", false);
  config.add("HrepMode", false);
  config.add("withFriction", true);

  // Create simple contact
  mc_rbdyn::Contact leftFootContact(const_cast<mc_rbdyn::Robots &>(*robotsList.get()), 0, 1, "LeftFoot", "AllGround");

  // polytope dimension depending on moments option
  config("withMoments") ? Rn::setDimension(6) : Rn::setDimension(3);

  // Testing part
  for(auto _ : state)
  {
    // Intentionally building job and input in the test loop to emulate lib runtime
    mc_dynamic_polytopes::ContactPolytopeJob contactJob;
    contactJob.contactRBDyn_ = &leftFootContact;
    contactJob.ctl_ = NULL;

    // "job" inputs (simple defaults)
    auto & input = contactJob.input();
    input.contactName = "LeftFoot";
    input.mb = robot.mb();
    input.mbc = robot.mbc();
    input.tl = robot.tl();
    input.tu = robot.tu();
    input.surface = robot.surface("LeftFoot").copy();
    input.surfacePose = input.surface->X_0_s(robot);
    input.accW = robot.accW();
    input.refContactTransform = sva::PTransformd::Identity();
    input.frictionCoefficient = 0.7;
    input.forceScalingFactor = 1.0;

    // run computation manually
    contactJob.computeJob();
  }
}
// Register the function as a benchmark
BENCHMARK(BM_PolytopeJVRC1_3D_6DOF)->Unit(benchmark::kMillisecond);

static void BM_PolytopeJVRC1_6D_6DOF(benchmark::State & state)
{
  // Setup part (not timed)
  auto robots = mc_rbdyn::loadRobot(*mc_rbdyn::RobotLoader::get_robot_module("JVRC1"));
  auto & robot = robots->robot();
  auto envs = mc_rbdyn::loadRobot(*mc_rbdyn::RobotLoader::get_robot_module("env/ground"));
  auto & ground = envs->robot();

  mc_rbdyn::RobotsPtr robotsList(mc_rbdyn::Robots::make());
  robotsList->robotCopy(robot, robot.name());
  robotsList->robotCopy(ground, ground.name());

  mc_rtc::Configuration config;
  // Add required configuration keys if needed
  config.add("possibleContacts", std::vector<std::string>{"LeftFoot", "RightFoot"});
  config.add("withMoments", true);
  config.add("computeRegions", false);
  config.add("HrepMode", false);
  config.add("withFriction", true);

  // Create simple contact
  mc_rbdyn::Contact leftFootContact(const_cast<mc_rbdyn::Robots &>(*robotsList.get()), 0, 1, "LeftFoot", "AllGround");

  // polytope dimension depending on moments option
  config("withMoments") ? Rn::setDimension(6) : Rn::setDimension(3);

  // Testing part
  for(auto _ : state)
  {
    // Intentionally building job and input in the test loop to emulate lib runtime
    mc_dynamic_polytopes::ContactPolytopeJob contactJob;
    contactJob.contactRBDyn_ = &leftFootContact;
    contactJob.ctl_ = NULL;

    // "job" inputs (simple defaults)
    auto & input = contactJob.input();
    input.contactName = "LeftFoot";
    input.mb = robot.mb();
    input.mbc = robot.mbc();
    input.tl = robot.tl();
    input.tu = robot.tu();
    input.surface = robot.surface("LeftFoot").copy();
    input.surfacePose = input.surface->X_0_s(robot);
    input.accW = robot.accW();
    input.refContactTransform = sva::PTransformd::Identity();
    input.frictionCoefficient = 0.7;
    input.forceScalingFactor = 1.0;

    // run computation manually
    contactJob.computeJob();
  }
}

BENCHMARK(BM_PolytopeJVRC1_6D_6DOF)->Unit(benchmark::kMillisecond);
