/* @todo
- regarder instanciation faux robot comme dans tests unitaires mc_rtc
- vérifier à quel point nb total de dofs impacte (taille calcul jacobienne inertie etc)
- déterminer runtime loop quand le future() du thread de la lib est récupéré
*/

#include <mc_rbdyn/RobotLoader.h>
#include <mc_rtc/Configuration.h>
#include <benchmark/benchmark.h>
#include <mc_dynamic_polytopes/DynamicPolytope.h>

static void BM_PolytopeJVRC1(benchmark::State & state)
{
  // Setup part (not timed)
  auto robots = mc_rbdyn::loadRobot(*mc_rbdyn::RobotLoader::get_robot_module("JVRC1"));
  auto & robot = robots->robot();
  auto envs = mc_rbdyn::loadRobot(*mc_rbdyn::RobotLoader::get_robot_module("env/ground"));
  auto & ground = envs->robot();

  mc_rbdyn::RobotsPtr robotsList(mc_rbdyn::Robots::make());
  robotsList->robotCopy(robot, robot.name());
  robotsList->robotCopy(ground, ground.name());
  // robotsList->load()

  mc_rtc::Configuration config;
  // Add required configuration keys if needed
  config.add("possibleContacts", std::vector<std::string>{"LeftFoot", "RightFoot"});
  config.add("withMoments", false);
  config.add("computeRegions", true);
  config.add("HrepMode", false);
  config.add("withFriction", true);
  mc_dynamic_polytopes::DynamicPolytope polytopeComputer("test", robot, config);

  // Set current active contacts (normally in controller run)
  std::map<std::string, mc_rbdyn::Contact &> Contacts; // const_cast<mc_rbdyn::Robots &>(*robotsList)
  mc_rbdyn::Contact leftFootContact(const_cast<mc_rbdyn::Robots &>(*robotsList.get()), 0, 1, "LeftFoot", "AllGround");
  Contacts.emplace(leftFootContact.r1Surface()->name(), leftFootContact);

  // Testing part
  for(auto _ : state)
  {
    polytopeComputer.setControllerContacts(Contacts);
    polytopeComputer.computeRegions();
    // @todo wait for polytope result
  }
}
// Register the function as a benchmark
BENCHMARK(BM_PolytopeJVRC1);
