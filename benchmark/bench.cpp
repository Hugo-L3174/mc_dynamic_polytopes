/* @todo
- Check if total number of dofs in the robot has an impact (size of jacobian and inertia matrices)
- Check if need to set random posture and dynamics to robot ? some calculations may be optimized if full zeros
*/

#include <mc_rbdyn/RobotLoader.h>
#include <mc_rtc/Configuration.h>
#include <benchmark/benchmark.h>
#include <mc_dynamic_polytopes/DynamicPolytope.h>

void addDummySurface(mc_rbdyn::Robot & robot, std::string surfaceName, std::string parentBodyName)
{
  // create a dummy square surface
  std::vector<std::pair<double, double>> points = {{-0.1, -0.1}, {0.1, -0.1}, {0.1, 0.1}, {-0.1, 0.1}};
  std::shared_ptr<mc_rbdyn::Surface> surfacePtr(
      new mc_rbdyn::PlanarSurface(surfaceName, parentBodyName, sva::PTransformd::Identity(), "plastic", points));
  robot.addSurface(surfacePtr);
}

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
BENCHMARK(BM_PolytopeJVRC1_3D_6DOF)->Unit(benchmark::kMillisecond);
// @todo consider setting number of repetitons to get standard deviation with
// ->Repetitions(30)->DisplayAggregatesOnly(true);

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

BENCHMARK(BM_PolytopeJVRC1_6D_6DOF)->Unit(benchmark::kMillisecond);

static void BM_PolytopeJVRC1_3D_7DOF(benchmark::State & state)
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
  // JVRC1 arm: WAIST_Y -> WAIST_P -> WAIST_R -> R_SHOULDER_P -> R_SHOULDER_R -> R_SHOULDER_Y -> R_ELBOW_P
  // ==> contact on R_ELBOW_P link is 7 dofs
  addDummySurface(robot, "RForearm", "R_ELBOW_P_S");
  mc_rbdyn::Contact RForearmContact(*robots, 0, 1, "RForearm", "AllGround");

  // polytope dimension depending on moments option
  config("withMoments") ? Rn::setDimension(6) : Rn::setDimension(3);

  // Testing part
  for(auto _ : state)
  {
    // Intentionally building job and input in the test loop to emulate lib runtime
    mc_dynamic_polytopes::ContactPolytopeJob contactJob;
    contactJob.contactRBDyn_ = &RForearmContact;
    contactJob.ctl_ = NULL;
    contactJob.HrepMode_ = config("HrepMode");
    contactJob.combineWithFriction_ = config("withFriction");

    // "job" inputs (simple defaults)
    auto & input = contactJob.input();
    input.contactName = "RForearm";
    input.mb = robot.mb();
    input.mbc = robot.mbc();
    input.tl = robot.tl();
    input.tu = robot.tu();
    input.surface = robot.surface("RForearm").copy();
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
BENCHMARK(BM_PolytopeJVRC1_3D_7DOF)->Unit(benchmark::kMillisecond);

static void BM_PolytopeJVRC1_6D_7DOF(benchmark::State & state)
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
  // JVRC1 arm: WAIST_Y -> WAIST_P -> WAIST_R -> R_SHOULDER_P -> R_SHOULDER_R -> R_SHOULDER_Y -> R_ELBOW_P
  // ==> contact on R_ELBOW_P link is 7 dofs
  addDummySurface(robot, "RForearm", "R_ELBOW_P_S");
  mc_rbdyn::Contact RForearmContact(*robots, 0, 1, "RForearm", "AllGround");

  // polytope dimension depending on moments option
  config("withMoments") ? Rn::setDimension(6) : Rn::setDimension(3);

  // Testing part
  for(auto _ : state)
  {
    // Intentionally building job and input in the test loop to emulate lib runtime
    mc_dynamic_polytopes::ContactPolytopeJob contactJob;
    contactJob.contactRBDyn_ = &RForearmContact;
    contactJob.ctl_ = NULL;
    contactJob.HrepMode_ = config("HrepMode");
    contactJob.combineWithFriction_ = config("withFriction");

    // "job" inputs (simple defaults)
    auto & input = contactJob.input();
    input.contactName = "RForearm";
    input.mb = robot.mb();
    input.mbc = robot.mbc();
    input.tl = robot.tl();
    input.tu = robot.tu();
    input.surface = robot.surface("RForearm").copy();
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
BENCHMARK(BM_PolytopeJVRC1_6D_7DOF)->Unit(benchmark::kMillisecond);

static void BM_PolytopeJVRC1_3D_8DOF(benchmark::State & state)
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
  // JVRC1 arm: WAIST_Y -> WAIST_P -> WAIST_R -> R_SHOULDER_P -> R_SHOULDER_R -> R_SHOULDER_Y -> R_ELBOW_P -> R_ELBOW_Y
  // ==> contact on R_ELBOW_Y link is 8 dofs
  addDummySurface(robot, "RWrist", "R_ELBOW_Y_S");
  mc_rbdyn::Contact RWristContact(*robots, 0, 1, "RWrist", "AllGround");

  // polytope dimension depending on moments option
  config("withMoments") ? Rn::setDimension(6) : Rn::setDimension(3);

  // Testing part
  for(auto _ : state)
  {
    // Intentionally building job and input in the test loop to emulate lib runtime
    mc_dynamic_polytopes::ContactPolytopeJob contactJob;
    contactJob.contactRBDyn_ = &RWristContact;
    contactJob.ctl_ = NULL;
    contactJob.HrepMode_ = config("HrepMode");
    contactJob.combineWithFriction_ = config("withFriction");

    // "job" inputs (simple defaults)
    auto & input = contactJob.input();
    input.contactName = "RWrist";
    input.mb = robot.mb();
    input.mbc = robot.mbc();
    input.tl = robot.tl();
    input.tu = robot.tu();
    input.surface = robot.surface("RWrist").copy();
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
BENCHMARK(BM_PolytopeJVRC1_3D_8DOF)->Unit(benchmark::kMillisecond);

static void BM_PolytopeJVRC1_6D_8DOF(benchmark::State & state)
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
  // JVRC1 arm: WAIST_Y -> WAIST_P -> WAIST_R -> R_SHOULDER_P -> R_SHOULDER_R -> R_SHOULDER_Y -> R_ELBOW_P -> R_ELBOW_Y
  // ==> contact on R_ELBOW_Y link is 8 dofs
  addDummySurface(robot, "RWrist", "R_ELBOW_Y_S");
  mc_rbdyn::Contact RWristContact(*robots, 0, 1, "RWrist", "AllGround");

  // polytope dimension depending on moments option
  config("withMoments") ? Rn::setDimension(6) : Rn::setDimension(3);

  // Testing part
  for(auto _ : state)
  {
    // Intentionally building job and input in the test loop to emulate lib runtime
    mc_dynamic_polytopes::ContactPolytopeJob contactJob;
    contactJob.contactRBDyn_ = &RWristContact;
    contactJob.ctl_ = NULL;
    contactJob.HrepMode_ = config("HrepMode");
    contactJob.combineWithFriction_ = config("withFriction");

    // "job" inputs (simple defaults)
    auto & input = contactJob.input();
    input.contactName = "RWrist";
    input.mb = robot.mb();
    input.mbc = robot.mbc();
    input.tl = robot.tl();
    input.tu = robot.tu();
    input.surface = robot.surface("RWrist").copy();
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
BENCHMARK(BM_PolytopeJVRC1_6D_8DOF)->Unit(benchmark::kMillisecond);

static void BM_PolytopeJVRC1_3D_9DOF(benchmark::State & state)
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
  // JVRC1 arm: WAIST_Y -> WAIST_P -> WAIST_R -> R_SHOULDER_P -> R_SHOULDER_R -> R_SHOULDER_Y -> R_ELBOW_P -> R_ELBOW_Y
  // -> R_WRIST_R
  // ==> contact on R_WRIST_R_S link is 9 dofs
  addDummySurface(robot, "RHand", "R_WRIST_R_S");
  mc_rbdyn::Contact RHandContact(*robots, 0, 1, "RHand", "AllGround");

  // polytope dimension depending on moments option
  config("withMoments") ? Rn::setDimension(6) : Rn::setDimension(3);

  // Testing part
  for(auto _ : state)
  {
    // Intentionally building job and input in the test loop to emulate lib runtime
    mc_dynamic_polytopes::ContactPolytopeJob contactJob;
    contactJob.contactRBDyn_ = &RHandContact;
    contactJob.ctl_ = NULL;
    contactJob.HrepMode_ = config("HrepMode");
    contactJob.combineWithFriction_ = config("withFriction");

    // "job" inputs (simple defaults)
    auto & input = contactJob.input();
    input.contactName = "RHand";
    input.mb = robot.mb();
    input.mbc = robot.mbc();
    input.tl = robot.tl();
    input.tu = robot.tu();
    input.surface = robot.surface("RHand").copy();
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
BENCHMARK(BM_PolytopeJVRC1_3D_9DOF)->Unit(benchmark::kMillisecond);

static void BM_PolytopeJVRC1_6D_9DOF(benchmark::State & state)
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
  // JVRC1 arm: WAIST_Y -> WAIST_P -> WAIST_R -> R_SHOULDER_P -> R_SHOULDER_R -> R_SHOULDER_Y -> R_ELBOW_P -> R_ELBOW_Y
  // -> R_WRIST_R
  // ==> contact on R_WRIST_R_S link is 9 dofs
  addDummySurface(robot, "RHand", "R_WRIST_R_S");
  mc_rbdyn::Contact RHandContact(*robots, 0, 1, "RHand", "AllGround");

  // polytope dimension depending on moments option
  config("withMoments") ? Rn::setDimension(6) : Rn::setDimension(3);

  // Testing part
  for(auto _ : state)
  {
    // Intentionally building job and input in the test loop to emulate lib runtime
    mc_dynamic_polytopes::ContactPolytopeJob contactJob;
    contactJob.contactRBDyn_ = &RHandContact;
    contactJob.ctl_ = NULL;
    contactJob.HrepMode_ = config("HrepMode");
    contactJob.combineWithFriction_ = config("withFriction");

    // "job" inputs (simple defaults)
    auto & input = contactJob.input();
    input.contactName = "RHand";
    input.mb = robot.mb();
    input.mbc = robot.mbc();
    input.tl = robot.tl();
    input.tu = robot.tu();
    input.surface = robot.surface("RHand").copy();
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
BENCHMARK(BM_PolytopeJVRC1_6D_9DOF)->Unit(benchmark::kMillisecond);
