#include "DynamicPolytope.h"
#include "SupportRegion.h"

DynamicPolytope::DynamicPolytope(const std::string & name,
                                 std::set<std::string> contactNames,
                                 const mc_rbdyn::Robot & robot)
: name_(fmt::format("DynamicPolytope_" + name)), possibleContacts_(contactNames),
  robotNetWrench_(sva::ForceVecd::Zero()), robot_(robot)
{
  // Init dimension
  Rn::setDimension(3);

  for(const auto contact : contactNames)
  {
    // init triangles maps
    std::vector<std::array<Eigen::Vector3d, 3UL>> newTrianglesArray;
    frictionConesTrianglesMap_.emplace(contact, newTrianglesArray);
    forcePolyTrianglesMap_.emplace(contact, newTrianglesArray);
    momentPolytopesTrianglesMap_.emplace(contact, newTrianglesArray);

    // init cones maps
    boost::shared_ptr<Polytope_Rn> newCone(new Polytope_Rn());
    frictionCones_.emplace(contact, newCone);

    boost::shared_ptr<Polytope_Rn> newPoly(new Polytope_Rn());
    forcePolytopes_.emplace(contact, newPoly);

    boost::shared_ptr<Polytope_Rn> newMomentCone(new Polytope_Rn());
    frictionConesMoments_.emplace(contact, newMomentCone);

    HRepX3d newPlanes;
    frictionConesPlanes_.emplace(contact, newPlanes);

    refContactPoses_.emplace(contact, robot_.surfacePose(contact));

    ContactTimers newTimers;
    contactsTimers_.emplace(contact, newTimers);
  }
  // init CWC polytope
  CWCForces_.reset(new Polytope_Rn());
  CWCMoments_.reset(new Polytope_Rn());

  // init zmp region and intersection with ecmp region
  zmpRegion_.reset(new Polytope_Rn());
  CWCMomentLessUnCstr_.reset(new Polytope_Rn());
  zeroMomentRegion_.reset(new Polytope_Rn());

  // Start a computing thread that will run continuously with the internal robot and contacts updated in the controller,
  // mutex protected
  mainComputeThread_ = std::thread(&DynamicPolytope::computeRegions, this);
#ifndef WIN32
  // Lower thread priority so that it has a lesser priority than the real time thread
  auto th_handle = mainComputeThread_.native_handle();
  int policy = 0;
  sched_param param{};
  pthread_getschedparam(th_handle, &policy, &param);
  param.sched_priority = 10;
  if(pthread_setschedparam(th_handle, SCHED_RR, &param) != 0)
  {
    // mc_rtc::log::warning("[{}] Failed to lower thread priority. If you are running on a real-time system, this might "
    //                      "cause latency to the real-time loop.",
    //                      name_);
  }
#endif
}

DynamicPolytope::~DynamicPolytope()
{
  stopThread();
}

void DynamicPolytope::load(const mc_rtc::Configuration & config)
{
  if(auto gui = config.find("polyhedronForce"))
  {
    polyForceConfig_.fromConfig(*gui);
  }
  if(auto gui = config.find("polyhedronMoment"))
  {
    polyMomentConfig_.fromConfig(*gui);
  }
  if(auto gui = config.find("polyhedronZMP"))
  {
    polyZMPConfig_.fromConfig(*gui);
  }
  if(auto gui = config.find("polyhedronZeroMomentArea"))
  {
    polyZeroMomentAreaConfig_.fromConfig(*gui);
  }
  config("withMoments", withMoments_);
}

void DynamicPolytope::addToLogger(mc_rtc::Logger & logger, const std::string & prefix)
{
  logger.addLogEntry("perf_" + prefix + "totalLoop", this, [this]() { return dt_loop_total().count(); });
  logger.addLogEntry("perf_" + prefix + "computeContactSet", this, [this]() { return dt_contactSet().count(); });
  logger.addLogEntry("perf_" + prefix + "minkSum", this, [this]() { return dt_minkSum().count(); });
  logger.addLogEntry("perf_" + prefix + "updatePlanes", this, [this]() { return dt_updatePlanes().count(); });
  logger.addLogEntry("perf_" + prefix + "guiTrianglesContacts", this,
                     [this]() { return dt_guiTrianglesContacts().count(); });
  logger.addLogEntry("perf_" + prefix + "guiTrianglesRegions", this,
                     [this]() { return dt_guiTrianglesRegions().count(); });
  logger.addLogEntry("perf_" + prefix + "zeroMomentIntersection", this,
                     [this]() { return dt_zeroMomentInter().count(); });
  for(const auto & contact : possibleContacts_)
  {
    logger.addLogEntry("perf_" + prefix + contact + "_frictionCone", this,
                       [this, contact]() { return getContact_dt_frictionCone(contact).count(); });
    logger.addLogEntry("perf_" + prefix + contact + "_forcePolytope", this,
                       [this, contact]() { return getContact_dt_forcePolytope(contact).count(); });
    logger.addLogEntry("perf_" + prefix + contact + "_intersection", this,
                       [this, contact]() { return getContact_dt_intersection(contact).count(); });
    logger.addLogEntry("perf_" + prefix + contact + "_Total", this,
                       [this, contact]() { return getContact_dt_Total(contact).count(); });
  }
}

void DynamicPolytope::computeRegions()
{
  // This section allows to block main thread until cv_ signals that computation can start, then loops conditionally on
  // the computing_ variable that will also be used for thread stopping.
  std::mutex mutex;
  std::unique_lock<std::mutex> lk(mutex);
  cv_.wait(lk, [this]() -> bool { return computing_; });

  while(computing_)
  {
    auto start_loop = mc_rtc::clock::now();
    // Lock contact set mutex, set active contacts to be used in internal loops then unlock mutex
    setCurrentContacts();

    // Step 1.1: launch individual force polytopes and friction cones calculations in separate threads as they are
    // independant, then their intersection
    computeFeasibleForcesFromContactSet(robot_);

    // Step 1.2: launch ZMP region calculation at the same time, also independant
    if(!zmpThread_.joinable())
    {
      zmpThread_ = std::thread(&DynamicPolytope::computeZMPRegion, this, std::ref(robot_), std::ref(contactsRBDyn_));
    }

    // Step 1.3: launch CWC Momentlesss calculation at the same time, also independant
    if(!cwcMomentLessThread_.joinable())
    {
      cwcMomentLessThread_ = std::thread(&DynamicPolytope::computeMomentLessForceCone, this, std::ref(robot_), std::ref(contactsRBDyn_));
    }

    // Step 2: wait for finished individual feasible regions
    // Wait for finished threads and join them, then update their matrix constraints
    for(const auto contactName : activeContacts_)
    {
      feasiblePolytopesThreadsMutex_.lock();
      feasiblePolytopesThreads_.at(contactName).join();
      feasiblePolytopesThreads_.erase(contactName);
      feasiblePolytopesThreadsMutex_.unlock();

      std::lock_guard<std::mutex> lock(getContactMutex(frictionConesPlanesMutexes_, contactName));
      updatePlanesMatrixConstraint(frictionCones_.at(contactName), frictionConesPlanes_.at(contactName).first,
                                   frictionConesPlanes_.at(contactName).second);
    }

    dt_compute_contactSet_ = mc_rtc::clock::now() - start_loop;

    // Update contacts GUI
    updateTrianglesContactsGUIPolitopix();

    // Steps 3-6: launch the rest everytime the previous full region was computed
    if(VRPRegionComputed_)
    {
      if(minkSumThread_.joinable())
      {
        minkSumThread_.join();
      }

      VRPRegionComputed_ = false;
      minkSumThread_ = std::thread(&DynamicPolytope::computeVRPRegionWithMinkSum, this);
    }

    // Regions GUI is now updated at the end of their thread to ensure scaling and translation are done before updating
    dt_loop_total_ = mc_rtc::clock::now() - start_loop;
  }
}

void DynamicPolytope::buildFrictionConeFromContactWithVrep(int numberOfFrictionSides,
                                                           const sva::PTransformd contactSurface,
                                                           boost::shared_ptr<Polytope_Rn> & frictionCone,
                                                           std::mutex & frictionConeMutex,
                                                           double m_frictionCoef,
                                                           double maxForce)
{
  int dim = 3;
  boost::shared_ptr<Polytope_Rn> newCone(new Polytope_Rn());
  // for now generate cone generates only the directions for the rays: we assume it is a polyhedral cone
  auto generators =
      generatePolyhedralConeGens(numberOfFrictionSides, contactSurface.rotation(), m_frictionCoef, maxForce);
  // here we manipulate polytope objects so need to add origin as a generator on the polyhedral cone
  generators.emplace_back(Eigen::Vector3d::Zero());
  for(const auto g : generators)
  {
    boost::shared_ptr<Generator_Rn> gn(new Generator_Rn(dim));
    boost::numeric::ublas::vector<double> coords(3);
    coords.insert_element(0, g.x());
    coords.insert_element(1, g.y());
    coords.insert_element(2, g.z());
    gn->setCoordinates(coords);
    newCone->addGenerator(gn);
    // mc_rtc::log::info("Creating cone with vertex {}", g.transpose());
  }
  // update faces of the cone
  // The face computations are necessary for the minkowsky sum using normal fans
  DoubleDescriptionFromGenerators::Compute(newCone, 1000);
  // lock cone mutex, then reset cone pointer to newly computed cone
  std::lock_guard<std::mutex> lock(frictionConeMutex);
  frictionCone.reset();
  frictionCone = newCone;
  // mc_rtc::log::info("Created cone of dim {} with {} generators", forceCone->dimension(),
  //                   forceCone->numberOfGenerators());
}

void DynamicPolytope::buildFrictionConeFromContactWithHrep(int numberOfFrictionSides,
                                                           const sva::PTransformd contactSurface,
                                                           boost::shared_ptr<Polytope_Rn> & frictionCone,
                                                           std::mutex & frictionConeMutex,
                                                           double m_frictionCoef)
{
  int dim = 3;
  boost::shared_ptr<Polytope_Rn> newCone(new Polytope_Rn());

  // get the friction cones planes
  auto Hrep = generatePolyhedralConeHRep(numberOfFrictionSides, contactSurface.rotation(), m_frictionCoef);
  // mc_rtc::log::info("Hrep dims are {} rows, {} columns", Hrep.rows(), Hrep.cols());

  // Adding the planes as the H-representation of the cone directly
  for(size_t i = 0; i < Hrep.rows(); i++)
  {
    boost::shared_ptr<HalfSpace_Rn> hs(new HalfSpace_Rn(dim));
    boost::numeric::ublas::vector<double> normal(3);
    normal.insert_element(0, Hrep.row(i).coeff(0));
    normal.insert_element(1, Hrep.row(i).coeff(1));
    normal.insert_element(2, Hrep.row(i).coeff(2));
    hs->setCoefficients(normal);
    hs->setConstant(0.0);
    newCone->addHalfSpace(hs);
  }

  // XXX No need for DD anymore because polyhedron is unbounded, it is a polyhedral cone
  // politopixAPI::computeDoubleDescriptionWithoutCheck(newCone, 10000);

  // Computing double description for polyhedral cone (polytope would cap the unbounded part)
  // TODO write DD function for polyhedralCones in politopix and switch friction cones type

  // mc_rtc::log::info("Computed friction cone with {} hs and {} gens", newCone->numberOfHalfSpaces(),
  //                   newCone->numberOfGenerators());
  // lock cone mutex, then reset cone pointer to newly computed cone
  std::lock_guard<std::mutex> lock(frictionConeMutex);
  frictionCone.reset();
  frictionCone = newCone;
}

void DynamicPolytope::buildWrenchConeFromContact(int numberOfFrictionSides,
                                                 std::pair<std::pair<double, double>, sva::PTransformd> contactSurface,
                                                 boost::shared_ptr<Polytope_Rn> & forceCone,
                                                 boost::shared_ptr<Polytope_Rn> & momentPoly,
                                                 double m_frictionCoef,
                                                 double maxForce,
                                                 Eigen::Vector3d CoM)
{
  int dim = 3;
  boost::shared_ptr<Polytope_Rn> newForcePoly(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> newMomentPoly(new Polytope_Rn());
  // newCone
  // for now generate cone generates only the directions for the rays: we assume it is a polyhedral cone
  auto generators =
      generatePolyhedralConeGens(numberOfFrictionSides, contactSurface.second.rotation(), m_frictionCoef, maxForce);
  // contactPoint first is the pair of xHalfLength and yHalfLength of the rectangular contact
  std::vector<Eigen::Vector3d> points;
  points.emplace_back(contactSurface.first.first, contactSurface.first.second, 0);
  points.emplace_back(contactSurface.first.first, -contactSurface.first.second, 0);
  points.emplace_back(-contactSurface.first.first, -contactSurface.first.second, 0);
  points.emplace_back(-contactSurface.first.first, contactSurface.first.second, 0);

  // here we manipulate polytope objects so need to add origin as a generator on the polyhedral cone
  generators.emplace_back(Eigen::Vector3d::Zero());
  for(auto g : generators)
  {
    // this is the translational part: no variation of force depending on application point
    Eigen::Vector3d newForce = g;

    boost::shared_ptr<Generator_Rn> forceGN(new Generator_Rn(dim));
    boost::numeric::ublas::vector<double> forceCoords(3);
    forceCoords.insert_element(0, newForce.x());
    forceCoords.insert_element(1, newForce.y());
    forceCoords.insert_element(2, newForce.z());
    forceGN->setCoordinates(forceCoords);
    newForcePoly->addGenerator(forceGN);
    for(auto p : points)
    {
      // r is the extremity point of the surface : center + offset using half dimensions (-application point in world
      // to get desired frame)
      Eigen::Vector3d r = contactSurface.second.translation() + p - CoM;

      // here compute the generators using the "limits" of the surface contact (points r) and associate each generator
      // found for the contact this is the angular part of the 6d vector (resulting moment of the force generator at
      // application point, here CoM)
      Eigen::Vector3d newMoment = skewMatrix(r) * g;

      boost::shared_ptr<Generator_Rn> momentGN(new Generator_Rn(dim));
      boost::numeric::ublas::vector<double> momentCoords(3);

      momentCoords.insert_element(0, newMoment.x());
      momentCoords.insert_element(1, newMoment.y());
      momentCoords.insert_element(2, newMoment.z());
      momentGN->setCoordinates(momentCoords);
      newMomentPoly->addGenerator(momentGN);
    }
  }

  DoubleDescriptionFromGenerators::Compute(newForcePoly, 1000);
  DoubleDescriptionFromGenerators::Compute(newMomentPoly, 1000);

  forceCone->reset();
  forceCone = newForcePoly;
  momentPoly->reset();
  momentPoly = newMomentPoly;
}

void DynamicPolytope::buildActuationPolytopeFromContact(const std::string contactName,
                                                        const mc_rbdyn::Robot & robot,
                                                        boost::shared_ptr<Polytope_Rn> & actuationPolytope,
                                                        std::mutex & forceConeMutex)
{
  int dim = 3;
  boost::shared_ptr<Polytope_Rn> newPoly(new Polytope_Rn());

  /* Questions to answer:
  - Do I need 6d then use only half, or can I have 3d from the beginning?
  - Do I need sparse matrix or can I use dense matrix? probably dense is ok since decoupled
    Compute with full then check if result is different with dense
  */
  const int n_var = 6;
  // Removing underactuated dofs
  const int jacSize = robot.mb().nrDof() - 6;

  // Correct jacobian building needs names because we need to link the surface to a body and compute the jac of the body
  // + transform between surface and the body
  const mc_rbdyn::Surface & surface = robot.surface(contactName);
  std::string bodyName = surface.bodyName();
  sva::PTransformd X_s_b;

  X_s_b = robot.bodyPosW(bodyName) * (robot.surfacePose(contactName).inv());

  // Building jacobian to the contact frame
  rbd::Jacobian jac(robot.mb(), bodyName, X_s_b.inv().translation());
  // Dense jacobian
  Eigen::MatrixXd denseJac = jac.bodyJacobian(robot.mb(), robot.mbc());
  // mc_rtc::log::info("Dense Jacobian: \n{}", denseJac);

  // rotate + transpose dense jacobian
  // Eigen::MatrixXd denseJacT = (sva::PTransformd(X_s_b.rotation()).matrix() * denseJac).transpose();
  // mc_rtc::log::info("Dense Jacobian transposed + rotated: \n{}", denseJacT);

  // Allocate then fill sparse jacobian
  Eigen::MatrixXd fullJac = Eigen::MatrixXd::Zero(6, robot.mb().nrDof());
  jac.fullJacobian(robot.mb(), denseJac, fullJac);
  // mc_rtc::log::info("Sparse Jacobian: \n{}", fullJac);

  // rotate to contact frame + transpose for correct calculation
  const Eigen::MatrixXd fullJacT = (sva::PTransformd(X_s_b.rotation()).matrix() * fullJac).transpose();
  // mc_rtc::log::info("Sparse jacobian transposed and rotated to surface frame: \n{}", fullJacT);

  // Computing dynamics terms of the equation of motion
  rbd::ForwardDynamics forwardDyn(robot.mb());
  // Inertia
  forwardDyn.computeH(robot.mb(), robot.mbc());
  Eigen::MatrixXd inertiaMat = forwardDyn.H();

  // Coriolis
  // XXX should the coriolis matrix be used or the C vector of the forward dynamics?
  rbd::Coriolis coriolis(robot.mb());
  Eigen::MatrixXd coriolisMat = coriolis.coriolis(robot.mb(), robot.mbc());

  // XXX this version contains gravity and external forces, to check
  forwardDyn.computeC(robot.mb(), robot.mbc());
  Eigen::VectorXd coriolisVec = forwardDyn.C();

  // Joint vectors
  const Eigen::VectorXd qdot = rbd::dofToVector(robot.mb(), robot.mbc().alpha);
  const Eigen::VectorXd qddot = rbd::dofToVector(robot.mb(), robot.mbc().alphaD);

  // Aineq matrix is jacobian transpose block to transform wrench vec to joint torque (x2 with second one negative
  // for being over the lower torque limits)
  // We also take only the bottom blocks without top 6 to remove underactuated part
  Eigen::MatrixXd Aineq(2 * jacSize, n_var);
  Aineq.block(0, 0, jacSize, n_var) = fullJacT.block(6, 0, jacSize, 6);
  Aineq.block(jacSize, 0, jacSize, n_var) = -fullJacT.block(6, 0, jacSize, 6);
  // mc_rtc::log::info("A ineq: \n{}", Aineq);

  // Torque limits vectors
  Eigen::VectorXd upperTorqueLims = rbd::dofToVector(robot.mb(), robot.tu()).segment(6, jacSize);
  Eigen::VectorXd lowerTorqueLims = rbd::dofToVector(robot.mb(), robot.tl()).segment(6, jacSize);

  // Gravity vector
  Eigen::VectorXd gravVector(jacSize);
  gravVector.fill(-9.81);

  // delta: might only need C and H from FD (already have all the elements), then remove top 6 rows
  // Inertia body transpose * body vel? + Inertia* qdotdot + coriolis+ g(q)
  // XXX check validity of this term
  Eigen::VectorXd delta = inertiaMat.bottomRows(jacSize) + coriolisVec.bottomRows(jacSize);
  // mc_rtc::log::info("delta:\n{}", delta.transpose());
  // mc_rtc::log::info("inertia is {} rows and {} cols, dofs are {}", delta.rows(), delta.cols(), robot.mb().nrDof());
  // mc_rtc::log::info("coriolis mat is {} rows and {} cols", coriolisMat.rows(), coriolisMat.cols());
  // mc_rtc::log::info("coriolis vec is {} rows", coriolisVec.size());

  // bineq vec is torque upper and lower limits (with negative lower limits) + delta (inertia, coriolis, gravity...)
  Eigen::VectorXd bineq(Aineq.rows());

  bineq.segment(0, jacSize) = delta + upperTorqueLims;
  bineq.segment(jacSize, jacSize) = (-1.0) * (delta + lowerTorqueLims);
  // mc_rtc::log::info("b ineq:\n{}", bineq.transpose());

  // Create a half space from every inequality
  for(size_t i = 0; i < jacSize * 2; i++)
  {
    // Add half space only if row is not null, ie only if this dof plays into the contact force (reduces nb of planes to
    // simplify in polytope)
    if(!Aineq.row(i).isZero())
    {
      boost::shared_ptr<HalfSpace_Rn> hs(new HalfSpace_Rn(dim));
      boost::numeric::ublas::vector<double> coefficients(3);
      // Setting coefficients as force elements of the ineq matrix (3 last columns)
      coefficients.insert_element(0, Aineq.coeff(i, 3));
      coefficients.insert_element(1, Aineq.coeff(i, 4));
      coefficients.insert_element(2, Aineq.coeff(i, 5));
      hs->setCoefficients(coefficients);
      hs->setConstant(bineq.coeff(i));
      newPoly->addHalfSpace(hs);
    }
  }
  auto start_DD = mc_rtc::clock::now();
  // Compute double description from half spaces (not generators -> truncation with bounding box)
  politopixAPI::computeDoubleDescriptionWithoutCheck(newPoly, 2000);
  mc_rtc::duration_ms end_DD = mc_rtc::clock::now() - start_DD;
  // mc_rtc::log::info("time to run force poly DD : {}ms", end_DD.count());

  // lock cone mutex, then reset cone pointer to newly computed cone
  std::lock_guard<std::mutex> lock(forceConeMutex);
  actuationPolytope.reset();
  actuationPolytope = newPoly;
}

void DynamicPolytope::buildFeasiblePolytopeFromContact(const std::string contactName,
                                                       const mc_rbdyn::Robot & robot,
                                                       const sva::PTransformd refContactPose,
                                                       int numberOfFrictionSides,
                                                       double frictionCoeff,
                                                       boost::shared_ptr<Polytope_Rn> & frictionCone,
                                                       std::mutex & frictionConeMutex,
                                                       boost::shared_ptr<Polytope_Rn> & actuationPolytope,
                                                       std::mutex & forcePolyMutex,
                                                       ContactTimers & timers)
{
  // sva::PTransformd contactPose = robot.surfacePose(contactName);
  auto start_forcePoly = mc_rtc::clock::now();

  // update the correct force polytope in the map
  // launching thread and emplacing it in the threads map
  forcePolyThreadsMutex_.lock();
  forcePolyThreads_.emplace(contactName,
                            std::thread(&DynamicPolytope::buildActuationPolytopeFromContact, this, contactName,
                                        std::ref(robot), std::ref(actuationPolytope), std::ref(forcePolyMutex)));
#ifndef WIN32
  // Lower thread priority so that it has a lesser priority than the real time thread
  auto th_handle = forcePolyThreads_.at(contactName).native_handle();
  int policy = 0;
  sched_param param{};
  pthread_getschedparam(th_handle, &policy, &param);
  param.sched_priority = 10;
  if(pthread_setschedparam(th_handle, SCHED_RR, &param) != 0)
  {
    // XXX Check if warning exists on real time kernel
    // mc_rtc::log::warning(
    //     "[{}] {} thread: failed to lower thread priority. If you are running on a real-time system, this might "
    //     "cause latency to the real-time loop.",
    //     name_, contactName);
  }
#endif
  forcePolyThreadsMutex_.unlock();

  auto start_frictionCone = mc_rtc::clock::now();
  // launching thread and emplacing it in the threads map
  frictionConesThreadsMutex_.lock();
  frictionConesThreads_.emplace(contactName, std::thread(&DynamicPolytope::buildFrictionConeFromContactWithHrep, this,
                                                         numberOfFrictionSides, refContactPose, std::ref(frictionCone),
                                                         std::ref(frictionConeMutex), frictionCoeff));
#ifndef WIN32
  // Lower thread priority so that it has a lesser priority than the real time thread
  th_handle = frictionConesThreads_.at(contactName).native_handle();
  pthread_getschedparam(th_handle, &policy, &param);
  param.sched_priority = 10;
  if(pthread_setschedparam(th_handle, SCHED_RR, &param) != 0)
  {
    // mc_rtc::log::warning(
    //     "[{}] {} thread: failed to lower thread priority. If you are running on a real-time system, this might "
    //     "cause latency to the real-time loop.",
    //     name_, contactName);
  }
#endif
  frictionConesThreadsMutex_.unlock();
  // Wait for end of computations
  frictionConesThreadsMutex_.lock();
  frictionConesThreads_.at(contactName).join();
  frictionConesThreads_.erase(contactName);
  frictionConesThreadsMutex_.unlock();

  timers.dt_frictionCone = mc_rtc::clock::now() - start_frictionCone;

  forcePolyThreadsMutex_.lock();
  forcePolyThreads_.at(contactName).join();
  forcePolyThreads_.erase(contactName);
  forcePolyThreadsMutex_.unlock();

  timers.dt_forcePolytope = mc_rtc::clock::now() - start_forcePoly;

  auto start_intersection = mc_rtc::clock::now();
  // Compute intersection
  std::lock_guard<mutex> lockFriction(getContactMutex(frictionConesMutexes_, contactName));
  std::lock_guard<mutex> lockForce(getContactMutex(forcePolyMutexes_, contactName));

  // intersect friction cone planes with force polytope into friction cone object
  // mc_rtc::log::info("Friction cone {} before intersection: {} hs and {} gens", contactName,
  //                   frictionCones_.at(contactName)->numberOfHalfSpaces(),
  //                   frictionCones_.at(contactName)->numberOfGenerators());
  // mc_rtc::log::info("Force poly {} before intersection: {} hs and {} gens", contactName,
  //                   forcePolytopes_.at(contactName)->numberOfHalfSpaces(),
  //                   forcePolytopes_.at(contactName)->numberOfGenerators());
  politopixAPI::computeIntersectionWithoutCheck(forcePolytopes_.at(contactName), frictionCones_.at(contactName));
  // mc_rtc::log::info("Friction cone {} after intersection: {} hs and {} gens", contactName,
  //                   frictionCones_.at(contactName)->numberOfHalfSpaces(),
  //                   frictionCones_.at(contactName)->numberOfGenerators());
  // mc_rtc::log::info("Force poly {} after intersection: {} hs and {} gens", contactName,
  //                   forcePolytopes_.at(contactName)->numberOfHalfSpaces(),
  //                   forcePolytopes_.at(contactName)->numberOfGenerators());
  timers.dt_intersection = mc_rtc::clock::now() - start_intersection;
  timers.dt_contactTotal = mc_rtc::clock::now() - start_forcePoly;
}

void DynamicPolytope::computeFrictionConesFromContactSet(const mc_rbdyn::Robot & robot)
{
  Rn::setDimension(3);
  auto frictionCoeff = 0.7;
  auto nbFrictionSides = 5;
  auto maxForce = 250.;
  auto CoM = robot.com();

  for(const auto contactName : activeContacts_)
  {
    sva::PTransformd contactPose = robot.surfacePose(contactName);
    if(!withMoments_)
    {
      // update the correct cone in the map
      // launching thread and emplacing it in the threads map
      frictionConesThreads_.emplace(
          contactName, std::thread(&DynamicPolytope::buildFrictionConeFromContactWithHrep, this, nbFrictionSides,
                                   contactPose, std::ref(frictionCones_.at(contactName)),
                                   std::ref(getContactMutex(frictionConesMutexes_, contactName)), frictionCoeff));
#ifndef WIN32
      // Lower thread priority so that it has a lesser priority than the real time thread
      auto th_handle = frictionConesThreads_.at(contactName).native_handle();
      int policy = 0;
      sched_param param{};
      pthread_getschedparam(th_handle, &policy, &param);
      param.sched_priority = 10;
      if(pthread_setschedparam(th_handle, SCHED_RR, &param) != 0)
      {
        // mc_rtc::log::warning(
        //     "[{}] {} thread: failed to lower thread priority. If you are running on a real-time system, this might "
        //     "cause latency to the real-time loop.",
        //     name_, contactName);
      }
#endif
    }
    else
    {
      // find limits of contact area for moment limits
      double newContactHalfLength;
      double newContactHalfWidth;
      findHalfWidthLength(robot.surface(contactName), newContactHalfWidth, newContactHalfLength);
      std::pair<std::pair<double, double>, sva::PTransformd> newContact(
          std::pair<double, double>(newContactHalfLength, newContactHalfWidth), contactPose);

      // TODO thread moments versions as well: add moment mutex + put mutexes as arguments
      buildWrenchConeFromContact(nbFrictionSides, newContact, frictionCones_.at(contactName),
                                 frictionConesMoments_.at(contactName), frictionCoeff, maxForce, CoM);
    }
  }
}

void DynamicPolytope::computeForcePolyFromContactSet(const mc_rbdyn::Robot & robot)
{

  for(const auto contactName : activeContacts_)
  {
    // update the correct force polytope in the map
    // launching thread and emplacing it in the threads map
    forcePolyThreads_.emplace(contactName,
                              std::thread(&DynamicPolytope::buildActuationPolytopeFromContact, this, contactName,
                                          std::ref(robot), std::ref(forcePolytopes_.at(contactName)),
                                          std::ref(getContactMutex(forcePolyMutexes_, contactName))));
#ifndef WIN32
    // Lower thread priority so that it has a lesser priority than the real time thread
    auto th_handle = forcePolyThreads_.at(contactName).native_handle();
    int policy = 0;
    sched_param param{};
    pthread_getschedparam(th_handle, &policy, &param);
    param.sched_priority = 10;
    if(pthread_setschedparam(th_handle, SCHED_RR, &param) != 0)
    {
      // XXX Check if warning exists on real time kernel
      // mc_rtc::log::warning(
      //     "[{}] {} thread: failed to lower thread priority. If you are running on a real-time system, this might "
      //     "cause latency to the real-time loop.",
      //     name_, contactName);
    }
#endif
  }
}

void DynamicPolytope::computeFeasibleForcesFromContactSet(const mc_rbdyn::Robot & robot)
{
  Rn::setDimension(3);
  auto frictionCoeff = 0.7;
  auto nbFrictionSides = 5;
  std::lock_guard<std::mutex> lock(contactSetMutex_);
  for(const auto & contactName : activeContacts_)
  {
    // launching contact computation
    std::lock_guard<std::mutex> lockThreadMap(feasiblePolytopesThreadsMutex_);
    feasiblePolytopesThreads_.emplace(
        contactName,
        std::thread(
            &DynamicPolytope::buildFeasiblePolytopeFromContact, this, contactName, std::ref(robot),
            refContactPoses_.at(contactName), nbFrictionSides, frictionCoeff, std::ref(frictionCones_.at(contactName)),
            std::ref(getContactMutex(frictionConesMutexes_, contactName)), std::ref(forcePolytopes_.at(contactName)),
            std::ref(getContactMutex(forcePolyMutexes_, contactName)), std::ref(contactsTimers_.at(contactName))));
#ifndef WIN32
    // Lower thread priority so that it has a lesser priority than the real time thread
    auto th_handle = feasiblePolytopesThreads_.at(contactName).native_handle();
    int policy = 0;
    sched_param param{};
    pthread_getschedparam(th_handle, &policy, &param);
    param.sched_priority = 10;
    if(pthread_setschedparam(th_handle, SCHED_RR, &param) != 0)
    {
      // XXX Check if warning exists on real time kernel
      // mc_rtc::log::warning(
      //     "[{}] {} thread: failed to lower thread priority. If you are running on a real-time system, this might "
      //     "cause latency to the real-time loop.",
      //     name_, contactName);
    }
#endif
  }
}

void DynamicPolytope::computeMinkowskySumPolitopix()
{
  auto start_minkSum = mc_rtc::clock::now();
  boost::shared_ptr<Polytope_Rn> newForcePoly(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> newMomentPoly(new Polytope_Rn());

  // putting it in vector form for library function
  std::vector<boost::shared_ptr<Polytope_Rn>> polytopesForces;
  std::vector<boost::shared_ptr<Polytope_Rn>> polytopesMoments;

  // protecting set of contact names
  contactSetMutex_.lock();
  for(const auto & active : activeContacts_)
  {
    boost::shared_ptr<Polytope_Rn> newContactPoly(new Polytope_Rn());
    // The friction + force intersection was overwritten in the force polytopes
    // Locking polytope mutex
    // Copying the contact polytopes to use for mink computation
    getContactMutex(forcePolyMutexes_, active).lock();
    politopixAPI::copyPolytope(forcePolytopes_.at(active), newContactPoly);
    getContactMutex(forcePolyMutexes_, active).unlock();

    polytopesForces.emplace_back(newContactPoly);
    if(withMoments_)
    {
      polytopesMoments.emplace_back(frictionConesMoments_.at(active));
    }
  }
  contactSetMutex_.unlock();

  try
  {
    MinkowskiSum Mink(polytopesForces, newForcePoly);
    std::lock_guard<std::mutex> lock(CWCMutex_);
    CWCForces_.reset();
    CWCForces_ = newForcePoly;
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
    CWCMoments_.reset();
    CWCMoments_ = newMomentPoly;
  }

  dt_compute_minkSum_ = mc_rtc::clock::now() - start_minkSum;
}

void DynamicPolytope::computeECMPRegion(Eigen::Vector3d comPosition, const mc_rbdyn::Robot & robot)
{
  /* We want to scale the force polytope according to the eCMP expression:
  eCMP = c - sumF/(m*(g/Dz)) --> eCMP = c - sumF*(Dz/mg)
  --> scale force polytope by - Deltaz/mg */

  // Several steps needed to mirror the polytope correctly (including the planes orientations):

  // First we scale with the POSITIVE scale because scaling only modifies the HS offsets, not their normals so the
  // polytope would be inverted but the planes inside out
  double scale = comPosition.z() / (robot.mass() * 9.81);
  bool ok = TopGeomTools::scalingFactor(CWCForces_, scale);
  // Then we negate the polytope to mirror the generators and the normals correctly
  CWCForces_->negate();

  // Finally translate it from origin to the robot CoM to get eCMP region
  boost::numeric::ublas::vector<double> CoM(3);
  CoM[0] = comPosition.x();
  CoM[1] = comPosition.y();
  CoM[2] = comPosition.z();
  TopGeomTools::translate(CWCForces_, CoM);
}

void DynamicPolytope::VRPtranslation(double deltaZ)
 {
  boost::numeric::ublas::vector<double> deltaZVector(3);
  deltaZVector[0] = 0.;
  deltaZVector[1] = 0.;
  deltaZVector[2] = deltaZ;
  std::lock_guard<std::mutex> CWCLock(CWCMutex_);
  std::lock_guard<std::mutex> zeroMomentLock(zeroMomentMutex_);
  std::lock_guard<std::mutex> zeroMomentUnconstLock(CWCMomentLessUnCstrTrianglesMutex_);
  TopGeomTools::translate(CWCForces_, deltaZVector);
  TopGeomTools::translate(zeroMomentRegion_, deltaZVector);
  TopGeomTools::translate(CWCMomentLessUnCstr_, deltaZVector);
}
 
void DynamicPolytope::computeZMPRegion(const mc_rbdyn::Robot & robot, const std::vector<mc_rbdyn::Contact> & contacts)
{

  const auto comPosition = robot.com();
  contactSetMutex_.lock();
  const auto static_hull = static_zmp_region(robot,contacts,Eigen::Vector3d::Zero(),Eigen::Vector3d{0,0,1})[0];
  contactSetMutex_.unlock();
  int dim = 3;
  std::lock_guard<std::mutex> lock(ZMPMutex_);
  zmpRegion_->reset();

  for(auto & pt : static_hull)
  {
    boost::shared_ptr<Generator_Rn> gn(new Generator_Rn(dim));
    boost::numeric::ublas::vector<double> coords(3);
    coords.insert_element(0, pt.x());
    coords.insert_element(1, pt.y());
    coords.insert_element(2, pt.z());
    gn->setCoordinates(coords);
    zmpRegion_->addGenerator(gn);
  }


  boost::shared_ptr<Generator_Rn> CoMgn(new Generator_Rn(dim));
  boost::numeric::ublas::vector<double> coords(3);
  coords.insert_element(0, comPosition.x());
  coords.insert_element(1, comPosition.y());
  coords.insert_element(2, comPosition.z());
  CoMgn->setCoordinates(coords);
  zmpRegion_->addGenerator(CoMgn);

  DoubleDescriptionFromGenerators::Compute(zmpRegion_, 1000);
}

void DynamicPolytope::computeMomentLessForceCone(const mc_rbdyn::Robot & robot, const std::vector<mc_rbdyn::Contact> & contacts)
{
  const Eigen::Vector3d comPosition = robot.com();
  contactSetMutex_.lock();
  const std::vector<Eigen::Vector3d> rays = momentless_force_cone(robot,contacts,comPosition,Eigen::Vector3d{0,0,1});
  contactSetMutex_.unlock();
  int dim = 3;
  CWCMomentLessUnCstr_->reset();

  const double scale = comPosition.z() / (robot.mass() * mc_rtc::constants::GRAVITY);
  // const double scale = 0.8;

  for(auto & r : rays)
  {
    const Eigen::Vector3d pt = comPosition - r * scale;
    boost::shared_ptr<Generator_Rn> gn(new Generator_Rn(dim));
    boost::numeric::ublas::vector<double> coords(3);
    coords.insert_element(0, pt.x());
    coords.insert_element(1, pt.y());
    coords.insert_element(2, pt.z());
    gn->setCoordinates(coords);
    CWCMomentLessUnCstr_->addGenerator(gn);
  }


  boost::shared_ptr<Generator_Rn> CoMgn(new Generator_Rn(dim));
  boost::numeric::ublas::vector<double> coords(3);
  coords.insert_element(0, comPosition.x());
  coords.insert_element(1, comPosition.y());
  coords.insert_element(2, comPosition.z());
  CoMgn->setCoordinates(coords);
  CWCMomentLessUnCstr_->addGenerator(CoMgn);

  DoubleDescriptionFromGenerators::Compute(CWCMomentLessUnCstr_, 1000);
}


void DynamicPolytope::computeZeroMomentIntersection()
{
  auto start_zeroMomentIntersection = mc_rtc::clock::now();
  zeroMomentMutex_.lock();
  zeroMomentRegion_->reset();
  // making a deep copy of the force polytope to use as base for intersection with zmp region
  // (avoids shared_ptr problems)
  CWCMutex_.lock();
  politopixAPI::copyPolytope(CWCForces_, zeroMomentRegion_);
  CWCMutex_.unlock();

  ZMPMutex_.lock();
  // politopixAPI::computeIntersection(CWCForces_, zmpRegion_, zeroMomentRegion_);
  politopixAPI::computeIntersectionWithoutCheck(zeroMomentRegion_, zmpRegion_);
  ZMPMutex_.unlock();
  zeroMomentMutex_.unlock();

  dt_zeroMoment_intersection_ = mc_rtc::clock::now() - start_zeroMomentIntersection;
}

void DynamicPolytope::computeVRPRegionWithMinkSum()
{
  // Step 3: All feasible regions computed, start minkowsky sum computation
  computeMinkowskySumPolitopix();

  // Step 4: wait for finished mink sum to convert to eCMP region (no thread needed)
  computeECMPRegion(robot_.com(), robot_);

  // Step 5: wait for finished ZMP region to start zero moment intersection with eCMP region
  zmpThread_.join();
  cwcMomentLessThread_.join();
  computeZeroMomentIntersection();

  // Step 6: translate eCMP region and zero-moment intersection to get VRP regions
  VRPtranslation(robot_.com().z());

  // Update VRP planes internal variables to be fetched by controller
  auto start_updatePlanes = mc_rtc::clock::now();
  VRPPlanesMutex_.lock();
  CWCMutex_.lock();
  updatePlanesMatrixConstraint(CWCForces_, DCMVRPPlanes_.first, DCMVRPPlanes_.second);
  CWCMutex_.unlock();
  VRPPlanesMutex_.unlock();

  // Update zero moment region planes to be fetched by controller
  zeroMomentPlanesMutex_.lock();
  zeroMomentMutex_.lock();
  updatePlanesMatrixConstraint(zeroMomentRegion_, zeroMomentPlanes_.first, zeroMomentPlanes_.second);
  zeroMomentMutex_.unlock();
  zeroMomentPlanesMutex_.unlock();
  dt_update_planes_ = mc_rtc::clock::now() - start_updatePlanes;

  // Update GUI display of regions
  updateTrianglesRegionsGUIPolitopix();
  VRPRegionComputed_ = true;
}

void DynamicPolytope::computeMomentsRegion(Eigen::Vector3d comPosition, const mc_rbdyn::Robot & robot)
{
  // We scale the moment polytope according to the expression of the difference between eCMP and ZMP:
  // eCMP = ZMP + 1/(m*(g+\ddot(c)_z)) * (tau_y, - tau_x, 0.)
  // scale moment polytope by 1/(m*(g+\ddot(c)_z))
  double scale = 1 / (robot.mass() * (9.81 + robot.comAcceleration().z()));
  bool ok = TopGeomTools::scalingFactor(CWCMoments_, scale);

  // TODO change coords from varignon (check) + check that corresponds to inside of eCMP region?
}

// Eigen::Vector3d projectPointInVRPRegion(Eigen::Vector3d testedPoint)
// {
//   constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(CWCForces_->getListOfGenerators());
//   for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
//     for (unsigned int j=0; j<iteGN.current()->numberOfFacets(); j++) {
//       boost::shared_ptr<HalfSpace_Rn> HS = iteGN.current()->getFacet(j);
//       boost::numeric::ublas::vector<double> projectedPoint;
//       double halfSpaceNorm = norm_2(HS->vect());
//       double disPoint2Hyp = HS->computeDistancePointHyperplane(iteGN.current()->vect(), projectedPoint,
//       halfSpaceNorm);
//       //std::cout << "d" << iteGN.currentIteratorNumber() << " = " << disPoint2Hyp << std::endl;
//       if (disPoint2Hyp > 0.25*TOL || disPoint2Hyp < -0.25*TOL)
//         isVeryClose = false;
//     }
//     if (isVeryClose == false) {
//       averagePoint /= iteGN.current()->numberOfFacets();
//       iteGN.current()->setCoordinates(averagePoint);
//     }
//   }
// }

// Eigen::Vector3d projectPointInPolytope(Eigen::Vector3d testedPoint, boost::shared_ptr<Polytope_Rn> & polytope)
// {
//   boost::numeric::ublas::vector<double> closestProjectedPoint;
//   constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(polytope->getListOfGenerators());
//   for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
//     for (unsigned int j=0; j<iteGN.current()->numberOfFacets(); j++) {
//       boost::shared_ptr<HalfSpace_Rn> HS = iteGN.current()->getFacet(j);
//       boost::numeric::ublas::vector<double> projectedPoint;
//       double halfSpaceNorm = norm_2(HS->vect());
//       double disPoint2Hyp = HS->computeDistancePointHyperplane(iteGN.current()->vect(), projectedPoint,
//       halfSpaceNorm);
//       //std::cout << "d" << iteGN.currentIteratorNumber() << " = " << disPoint2Hyp << std::endl;
//       if (disPoint2Hyp > 0.25*TOL || disPoint2Hyp < -0.25*TOL)
//         isVeryClose = false;
//     }
//     if (isVeryClose == false) {
//       averagePoint /= iteGN.current()->numberOfFacets();
//       iteGN.current()->setCoordinates(averagePoint);
//     }
//   }
//   return
// }

bool DynamicPolytope::checkGravityCenterInPolytope(boost::shared_ptr<Polytope_Rn> & polytope)
{
  boost::numeric::ublas::vector<double> gravCenter(3);
  TopGeomTools::gravityCenter(polytope, gravCenter);
  Point_Rn testPoint(gravCenter(0), gravCenter(1), gravCenter(2));

  int resultPolito = polytope->checkPoint(testPoint);
  bool retResultPolitopix = false;
  if(resultPolito == 1)
  {
    retResultPolitopix = true;
  }

  Eigen::MatrixX3d Normals;
  Eigen::VectorXd Offsets;
  updatePlanesMatrixConstraint(polytope, Normals, Offsets);
  Eigen::Vector3d testEigen(gravCenter(0), gravCenter(1), gravCenter(2));
  Eigen::VectorXd test = Normals * testEigen - Offsets;
  bool retResultEigen = true;
  for(int coeff = 0; coeff < test.size(); coeff++)
  {
    if(test(coeff) > 0.0)
    {
      retResultEigen = false;
    }
  }

  return retResultPolitopix && retResultEigen;
}

void DynamicPolytope::updateTrianglesContactsGUIPolitopix()
{
  auto start_guiTriangles = mc_rtc::clock::now();
  // Protecting set of names (threaded) then copying
  contactSetMutex_.lock();
  auto activeContactSet = activeContacts_;
  auto contactsToBeRemoved = contactsToRemove_;
  contactSetMutex_.unlock();

  for(const auto & contact : activeContactSet)
  {
    auto contactPose = robot_.surfacePose(contact).translation();
    getContactMutex(frictionConeTrianglesMutexes_, contact).lock();
    getContactMutex(frictionConesMutexes_, contact).lock();
    update3DPolyTrianglesPolitopix(frictionCones_.at(contact), frictionConesTrianglesMap_.at(contact), guiScale_,
                                   contactPose);
    getContactMutex(frictionConesMutexes_, contact).unlock();
    getContactMutex(frictionConeTrianglesMutexes_, contact).unlock();
    if(withMoments_)
    {
      getContactMutex(momentTrianglesMutexes_, contact).lock();
      update3DPolyTrianglesPolitopix(frictionConesMoments_.at(contact), momentPolytopesTrianglesMap_.at(contact),
                                     guiScale_, contactPose);
      getContactMutex(momentTrianglesMutexes_, contact).unlock();
    }
    getContactMutex(forcePolyTrianglesMutexes_, contact).lock();
    getContactMutex(forcePolyMutexes_, contact).lock();
    update3DPolyTrianglesPolitopix(forcePolytopes_.at(contact), forcePolyTrianglesMap_.at(contact), guiScale_,
                                   contactPose);
    getContactMutex(forcePolyMutexes_, contact).unlock();
    getContactMutex(forcePolyTrianglesMutexes_, contact).unlock();
  }
  for(const auto & contact : contactsToBeRemoved)
  {
    getContactMutex(frictionConeTrianglesMutexes_, contact).lock();
    frictionConesTrianglesMap_.at(contact).clear();
    getContactMutex(frictionConeTrianglesMutexes_, contact).unlock();
    if(withMoments_)
    {
      getContactMutex(momentTrianglesMutexes_, contact).lock();
      momentPolytopesTrianglesMap_.at(contact).clear();
      getContactMutex(momentTrianglesMutexes_, contact).unlock();
    }
    getContactMutex(forcePolyTrianglesMutexes_, contact).lock();
    forcePolyTrianglesMap_.at(contact).clear();
    getContactMutex(forcePolyTrianglesMutexes_, contact).unlock();
  }
  dt_compute_guiTrianglesContacts_ = mc_rtc::clock::now() - start_guiTriangles;
}

void DynamicPolytope::updateTrianglesRegionsGUIPolitopix()
{
  auto start_guiTriangles = mc_rtc::clock::now();
  contactSetMutex_.lock();
  auto noActiveContacts = activeContacts_.empty();
  contactSetMutex_.unlock();

  if(!noActiveContacts)
  {
    // gui scale for CWC should be 1, it is position space and not force space (because eCMP)
    // update6DPolyTrianglesPolitopix(CWC_, CWCMomentTriangles_, CWCForceTriangles_, guiScale_);
    CWCForceTrianglesMutex_.lock();
    CWCMutex_.lock();
    update3DPolyTrianglesPolitopix(CWCForces_, CWCForceTriangles_, 1);
    CWCMutex_.unlock();
    CWCForceTrianglesMutex_.unlock();

    // scale 1 here: already position space
    ZMPTrianglesMutex_.lock();
    ZMPMutex_.lock();
    update3DPolyTrianglesPolitopix(zmpRegion_, ZMPTriangles_, 1);
    ZMPMutex_.unlock();
    ZMPTrianglesMutex_.unlock();

    CWCMomentLessUnCstrTrianglesMutex_.lock();
    update3DPolyTrianglesPolitopix(CWCMomentLessUnCstr_, CWCMomentLessUnCstrTriangles_, 1);
    CWCMomentLessUnCstrTrianglesMutex_.unlock();

    zeroMomentTrianglesMutex_.lock();
    zeroMomentMutex_.lock();
    update3DPolyTrianglesPolitopix(zeroMomentRegion_, zeroMomentTriangles_, 1);
    zeroMomentMutex_.unlock();
    zeroMomentTrianglesMutex_.unlock();

    if(withMoments_)
    {
      CWCMomentTrianglesMutex_.lock();
      update3DPolyTrianglesPolitopix(CWCMoments_, CWCMomentTriangles_, guiScale_);
      CWCMomentTrianglesMutex_.unlock();
    }
  }
  else
  {
    CWCForceTrianglesMutex_.lock();
    CWCForceTriangles_.clear();
    CWCForceTrianglesMutex_.unlock();

    ZMPTrianglesMutex_.lock();
    ZMPTriangles_.clear();
    ZMPTrianglesMutex_.unlock();

    zeroMomentTrianglesMutex_.lock();
    zeroMomentTriangles_.clear();
    zeroMomentTrianglesMutex_.unlock();

    CWCMomentLessUnCstrTrianglesMutex_.lock();
    CWCMomentLessUnCstrTriangles_.clear();
    CWCMomentLessUnCstrTrianglesMutex_.unlock();

    if(withMoments_)
    {
      CWCMomentTrianglesMutex_.lock();
      CWCMomentTriangles_.clear();
      CWCMomentTrianglesMutex_.unlock();
    }
  }
  dt_compute_guiTrianglesRegions_ = mc_rtc::clock::now() - start_guiTriangles;
}

void DynamicPolytope::updatePlanesMatrixConstraint(const boost::shared_ptr<Polytope_Rn> & polytope,
                                                   Eigen::MatrixX3d & Normals,
                                                   Eigen::VectorXd & Offsets)
{
  // XXX might be interesting to use a vector of mc_rbdyn::Plane in the future ? but matrixXd would be easier to put in
  // blocks in a handmade QP where if planes we would need a for loop on every plane to fill the QP matrix
  auto nbOfPlanes = polytope->numberOfHalfSpaces();
  Normals.resize(nbOfPlanes, 3);
  Offsets.resize(nbOfPlanes);

  for(int halfSpaceIndex = 0; halfSpaceIndex < nbOfPlanes; halfSpaceIndex++)
  {
    const auto halfSpace = polytope->getHalfSpace(halfSpaceIndex);
    // pushing back inverted normals
    Normals.row(halfSpaceIndex) << -halfSpace->getCoefficient(0), -halfSpace->getCoefficient(1),
        -halfSpace->getCoefficient(2);
    Offsets(halfSpaceIndex) = halfSpace->getConstant();
  }
}

void DynamicPolytope::addToGUI(mc_rtc::gui::StateBuilder & gui, double guiScale, std::vector<std::string> category)
{
  guiScale_ = guiScale;
  category.push_back(name_);
  auto conesCat = category;
  conesCat.push_back("Friction cones");
  auto CWCCat = category;
  CWCCat.push_back("Contact Wrench Cone");

  for(const auto contact : possibleContacts_)
  {
    gui.addElement(this, conesCat,
                   mc_rtc::gui::Polyhedron(fmt::format(contact + " frictions"), polyForceConfig_,
                                           [this, contact]() { return getFrictionConesTriangles(contact); }),
                   mc_rtc::gui::Polyhedron(fmt::format(contact + " forces"), polyMomentConfig_,
                                           [this, contact]() { return getForcePolyTriangles(contact); }),
                   mc_rtc::gui::Polyhedron(fmt::format(contact + " moments"), polyMomentConfig_,
                                           [this, contact]() { return getContactMomentTriangles(contact); }));
  }

  gui.addElement(
      this, CWCCat,
      mc_rtc::gui::Polyhedron("CWC forces", polyForceConfig_, [this]() { return getCWCForceTriangles(); }),
      mc_rtc::gui::Polyhedron("CWC moments", polyMomentConfig_, [this]() { return getCWCMomentTriangles(); }),
      mc_rtc::gui::Polyhedron("CWC moments uncstr", polyMomentConfig_, [this]() { return getCWCMomentLess(); }),
      mc_rtc::gui::Polyhedron("ZMP area", polyZMPConfig_, [this]() { return getZMPTriangles(); }),
      mc_rtc::gui::Polyhedron("Zero moment region", polyZeroMomentAreaConfig_,
                              [this]() { return getZeroMomentTriangles(); }));

  mc_rtc::gui::ArrowConfig Arrow;
  Arrow.scale = guiScale_;
  gui.addElement(this, category,
                 mc_rtc::gui::Point3D("eCMP", mc_rtc::gui::PointConfig(mc_rtc::gui::Color{1.0, 0.0, 0.0}, 0.03),
                                      [this]() -> const Eigen::Vector3d & { return eCMP_; }),
                 mc_rtc::gui::Arrow(
                     "Moments", Arrow, [this]() -> const Eigen::Vector3d { return Eigen::Vector3d::Zero(); },
                     [this]() -> const Eigen::Vector3d { return robotNetWrench_.couple(); }));
}
