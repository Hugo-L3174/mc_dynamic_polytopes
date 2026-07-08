/* @todo
- Check if total number of dofs in the robot has an impact (size of jacobian and inertia matrices)
- Check if need to set random posture and dynamics to robot ? some calculations may be optimized if full zeros
*/

#include <mc_rbdyn/RobotLoader.h>
#include <mc_rtc/Configuration.h>
#include <benchmark/benchmark.h>
#include <mc_dynamic_polytopes/DynamicPolytope.h>

static void BM_PolytopeJVRC1_3D_6DOF(benchmark::State & state)
{
  // Setup part (not timed)
  auto robotModule = mc_rbdyn::RobotLoader::get_robot_module("JVRC1");
  auto ground = mc_rbdyn::RobotLoader::get_robot_module("env/ground");
  auto robots = mc_rbdyn::loadRobots({robotModule, ground});
  auto & robot = robots->robot(robotModule->name);

  mc_rtc::Configuration config;
  config.add("withMoments", false);
  config.add("HrepMode", false);
  config.add("withFriction", true);

  // Create simple contact
  // JVRC1 leg: L_HIP_P -> L_HIP_R -> L_HIP_Y -> L_KNEE -> L_ANKLE_R -> L_ANKLE_P
  // ==> 6 dofs
  mc_rbdyn::Contact leftFootContact(*robots, 0, 1, "LeftFoot", "AllGround");

  // polytope dimension depending on moments option
  config("withMoments") ? Rn::setDimension(6) : Rn::setDimension(3);

  // Testing part
  for(auto _ : state)
  {
    // Intentionally building job and input in the test loop to emulate lib runtime
    mc_dynamic_polytopes::ContactPolytopeJob contactJob;
    contactJob.contactRBDyn_ = &leftFootContact;
    contactJob.ctl_ = NULL;
    contactJob.HrepMode_ = config("HrepMode");
    contactJob.combineWithFriction_ = config("withFriction");

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
BENCHMARK(BM_PolytopeJVRC1_3D_6DOF)->Unit(benchmark::kMillisecond)->Repetitions(100)->DisplayAggregatesOnly(true);

static void BM_PolytopeJVRC1_6D_6DOF(benchmark::State & state)
{
  // Setup part (not timed)
  auto robotModule = mc_rbdyn::RobotLoader::get_robot_module("JVRC1");
  auto ground = mc_rbdyn::RobotLoader::get_robot_module("env/ground");
  auto robots = mc_rbdyn::loadRobots({robotModule, ground});
  auto & robot = robots->robot(robotModule->name);

  mc_rtc::Configuration config;
  config.add("withMoments", true);
  config.add("HrepMode", false);
  config.add("withFriction", true);

  // Create simple contact
  // JVRC1 leg: L_HIP_P -> L_HIP_R -> L_HIP_Y -> L_KNEE -> L_ANKLE_R -> L_ANKLE_P
  // ==> 6 dofs
  mc_rbdyn::Contact leftFootContact(*robots, 0, 1, "LeftFoot", "AllGround");

  // polytope dimension depending on moments option
  config("withMoments") ? Rn::setDimension(6) : Rn::setDimension(3);

  // Testing part
  for(auto _ : state)
  {
    // Intentionally building job and input in the test loop to emulate lib runtime
    mc_dynamic_polytopes::ContactPolytopeJob contactJob;
    contactJob.contactRBDyn_ = &leftFootContact;
    contactJob.ctl_ = NULL;
    contactJob.HrepMode_ = config("HrepMode");
    contactJob.combineWithFriction_ = config("withFriction");

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

BENCHMARK(BM_PolytopeJVRC1_6D_6DOF)->Unit(benchmark::kMillisecond)->Repetitions(100)->DisplayAggregatesOnly(true);
