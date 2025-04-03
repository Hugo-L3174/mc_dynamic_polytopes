#include "DynamicPolytope.h"

DynamicPolytope::DynamicPolytope(const std::string & name,
                                 const mc_rbdyn::Robot & robot,
                                 const mc_rtc::Configuration & dynamicPolyConfig)
: name_(fmt::format("DynamicPolytope_" + name)), robot_(robot), config_(dynamicPolyConfig)
{
  // Init dimension
  Rn::setDimension(3);
  Rn::setTolerance(1.e-07);

  possibleContacts_ = dynamicPolyConfig("possibleContacts", std::set<std::string>{"LeftFoot", "RightFoot"});
  withMoments_ = dynamicPolyConfig("withMoments", false);
  computeRegions_ = dynamicPolyConfig("computeRegions", true);
  HrepMode_ = dynamicPolyConfig("HrepMode", false);

  double defaultForceScale = 1;
  int defaultFrictionSides = 5;
  double defaultFrictionCoeff = 0.5;

  for(const auto contact : possibleContacts_)
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
    // Here initialize planes as simple friction cones to have a sane constraint at first iteration
    // The rotX_r1_r2 matrix would be the identity in a default case
    newPlanes.first =
        generatePolyhedralConeHRep(defaultFrictionSides, Eigen::Matrix3d::Identity(), defaultFrictionCoeff);
    newPlanes.second = Eigen::VectorXd::Zero(newPlanes.first.rows());
    frictionConesPlanes_.emplace(contact, newPlanes);
    forcePolyPlanes_.emplace(contact, newPlanes);

    refContactTransforms_.emplace(contact, sva::PTransformd::Identity());

    ContactTimers newTimers;
    newTimers.dt_contactTotal = mc_rtc::duration_ms::zero();
    newTimers.dt_forcePolytope = mc_rtc::duration_ms::zero();
    newTimers.dt_frictionCone = mc_rtc::duration_ms::zero();
    newTimers.dt_intersection = mc_rtc::duration_ms::zero();
    contactsTimers_.emplace(contact, newTimers);

    forceScalingFactors_.emplace(contact, defaultForceScale);
    frictionCoefficients_.emplace(contact, defaultFrictionCoeff);
    numbersOfFrictionSides_.emplace(contact, defaultFrictionSides);
  }
  // init CWC polytope
  CWCForces_.reset(new Polytope_Rn());
  CWCMoments_.reset(new Polytope_Rn());

  // init zmp region and intersection with ecmp region
  zmpRegion_.reset(new Polytope_Rn());
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
    mc_rtc::log::warning("[{}] Failed to lower thread priority. If you are running on a real-time system, this might "
                         "cause latency to the real-time loop.",
                         name_);
  }
#endif
}

DynamicPolytope::~DynamicPolytope()
{
  stopThread();
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
    if(computeRegions_)
    {
      if(!zmpThread_.joinable())
      {
        zmpThread_ = std::thread(&DynamicPolytope::computeZMPRegion, this, robot_.com());
      }
    }

    // Step 2: wait for finished individual feasible regions
    // Wait for finished threads and join them, then update their matrix constraints
    for(const auto contactName : activeContacts_)
    {
      feasiblePolytopesThreadsMutex_.lock();
      feasiblePolytopesThreads_.at(contactName).join();
      feasiblePolytopesThreads_.erase(contactName);
      feasiblePolytopesThreadsMutex_.unlock();

      std::lock_guard<std::mutex> lockFriction(getContactMutex(frictionConesPlanesMutexes_, contactName));
      updatePlanesMatrixConstraint(frictionCones_.at(contactName), frictionConesPlanes_.at(contactName).first,
                                   frictionConesPlanes_.at(contactName).second);
      std::lock_guard<std::mutex> lockFeasible(getContactMutex(forcePolyPlanesMutexes_, contactName));
      updatePlanesMatrixConstraint(forcePolytopes_.at(contactName), forcePolyPlanes_.at(contactName).first,
                                   forcePolyPlanes_.at(contactName).second);
      std::lock_guard<std::mutex> lockContactSet(contactSetMutex_); // lock contact set for rbdyn contacts accessor
      updateRBDynPolytopes(forcePolyPlanes_.at(contactName).first, forcePolyPlanes_.at(contactName).second,
                           contactsRBDyn_.at(contactName));
    }

    dt_compute_contactSet_ = mc_rtc::clock::now() - start_loop;

    // Update contacts GUI
    updateTrianglesContactsGUIPolitopix();

    // Steps 3-6: launch the rest everytime the previous full region was computed
    if(computeRegions_)
    {
      if(VRPRegionComputed_)
      {
        if(minkSumThread_.joinable())
        {
          minkSumThread_.join();
        }

        VRPRegionComputed_ = false;
        minkSumThread_ = std::thread(&DynamicPolytope::computeVRPRegionWithMinkSum, this);
      }
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
                                                           const sva::PTransformd X_r1_r2,
                                                           boost::shared_ptr<Polytope_Rn> & frictionCone,
                                                           std::mutex & frictionConeMutex,
                                                           double m_frictionCoef)
{
  int dim = 3;
  boost::shared_ptr<Polytope_Rn> newCone(new Polytope_Rn());

  // get the friction cones planes
  auto Hrep = generatePolyhedralConeHRep(numberOfFrictionSides, X_r1_r2.rotation(), m_frictionCoef);
  // mc_rtc::log::info("Hrep dims are {} rows, {} columns", Hrep.rows(), Hrep.cols());

  // Adding the planes as the H-representation of the cone directly
  // XXX POLITOPIX: since politopix convention seems to be inverted for Hrep inequalities
  // i.e. they check that a0 + a1*x1 + a2*x2 ... >= 0 to belong inside
  // we push inverted normals
  for(auto i = 0; i < Hrep.rows(); i++)
  {
    boost::shared_ptr<HalfSpace_Rn> hs(new HalfSpace_Rn(dim));
    boost::numeric::ublas::vector<double> normal(3);
    normal.insert_element(0, Hrep.row(i).coeff(0));
    normal.insert_element(1, Hrep.row(i).coeff(1));
    normal.insert_element(2, Hrep.row(i).coeff(2));
    hs->setCoefficients(-normal);
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
                                                        std::mutex & forcePolyMutex,
                                                        double forceScalingFactor)
{
  int dim = 3;
  boost::shared_ptr<Polytope_Rn> newPoly(new Polytope_Rn());

  const int n_var = 6;
  // Removing underactuated dofs
  const int jacSize = robot.mb().nrDof() - 6;

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

  // Correct jacobian building needs names because we need to link the surface to a body and compute the jac of the body
  // + transform between surface and the body
  const mc_rbdyn::Surface & surface = robot.surface(contactName);
  std::string bodyName = surface.bodyName();
  sva::PTransformd X_s_b = robot.bodyPosW(bodyName) * (robot.surfacePose(contactName).inv());

  // XXX check if need transformation between body and frame OR between 0 and frame !
  // XXX here I used the bodyJacobian to be already in the body frame but is this correct?
  // Building jacobian to the contact frame
  rbd::Jacobian jac(robot.mb(), bodyName, X_s_b.inv().translation());

  // Dense jacobian
  Eigen::MatrixXd denseJac = jac.bodyJacobian(robot.mb(), robot.mbc());
  Eigen::MatrixXd denseJacDot = jac.bodyJacobianDot(robot.mb(), robot.mbc());
  // mc_rtc::log::info("Dense Jacobian: \n{}", denseJac);
  // mc_rtc::log::info("Dense Jacobian {} is of dimensions {}", contactName, denseJac.rows(), denseJac.cols());

  // Allocate then fill sparse jacobian
  Eigen::MatrixXd globalFullJac = Eigen::MatrixXd::Zero(6, robot.mb().nrDof());
  jac.fullJacobian(robot.mb(), denseJac, globalFullJac);
  Eigen::MatrixXd globalFullJacDot = Eigen::MatrixXd::Zero(6, robot.mb().nrDof());
  jac.fullJacobian(robot.mb(), denseJacDot, globalFullJacDot);
  // mc_rtc::log::info("Sparse Jacobian: \n{}", globalFullJac);

  // rotate to contact frame
  const Eigen::MatrixXd contactFullJac = sva::PTransformd(X_s_b.rotation()).matrix() * globalFullJac;
  // XXX to check ! simple rotation if jacDot matrix?
  const Eigen::MatrixXd contactFullJacDot = sva::PTransformd(X_s_b.rotation()).matrix() * globalFullJacDot;
  // transpose for correct calculation
  const Eigen::MatrixXd contactFullJacT = contactFullJac.transpose();
  // mc_rtc::log::info("Sparse jacobian transposed and rotated to surface frame: \n{}", contactFullJacT);

  // Torque limits vectors
  Eigen::VectorXd upperTorqueLims = rbd::dofToVector(robot.mb(), robot.tu());
  Eigen::VectorXd lowerTorqueLims = rbd::dofToVector(robot.mb(), robot.tl());
  Eigen::VectorXd upperTorqueLimsNoFloatingBase = rbd::dofToVector(robot.mb(), robot.tu()).segment(6, jacSize);
  Eigen::VectorXd lowerTorqueLimsNoFloatingBase = rbd::dofToVector(robot.mb(), robot.tl()).segment(6, jacSize);

  // Joint vectors
  const Eigen::VectorXd qdot = rbd::dofToVector(robot.mb(), robot.mbc().alpha);
  const Eigen::VectorXd qddot = rbd::dofToVector(robot.mb(), robot.mbc().alphaD);

  // XXX tentative selection vector for contact by using the contact jacobian
  Eigen::VectorXd actuatorSelectionVector = Eigen::VectorXd::Zero(robot.mb().nrDof());
  int numberOfActuatorsPlaying = 0;

  for(int i = 0; i < robot.mb().nrDof(); i++)
  {
    if(!contactFullJac.col(i).isZero())
    {
      actuatorSelectionVector(i) = Eigen::Index(1);
      numberOfActuatorsPlaying++;
    }
  }
  // mc_rtc::log::info("Actuator vector for {} is \n{}", contactName, actuatorSelectionVector);
  // mc_rtc::log::info("Actuator seletion vector for {} is of size {}", contactName, numberOfActuatorsPlaying);

  /***************************************************/
  // Hrep version

  if(HrepMode_)
  {
    // Aineq matrix is jacobian transpose block to transform wrench vec to joint torque (x2 with second one negative
    // for being over the lower torque limits)
    // We also take only the bottom blocks without top 6 to remove underactuated part
    Eigen::MatrixXd Aineq(2 * jacSize, n_var);
    Aineq.block(0, 0, jacSize, n_var) = contactFullJacT.block(6, 0, jacSize, 6);
    Aineq.block(jacSize, 0, jacSize, n_var) = -contactFullJacT.block(6, 0, jacSize, 6);
    // mc_rtc::log::info("A ineq: \n{}", Aineq);

    // Inertia floating base * floating base acc + Inertia without coupling terms * qdotdot + (coriolis+ g(q))

    Eigen::VectorXd delta = inertiaMat.bottomRows(jacSize).leftCols(6) * robot.accW().vector()
                            + inertiaMat.bottomRows(jacSize).rightCols(jacSize) * qddot.segment(6, jacSize)
                            + coriolisVec.bottomRows(jacSize);
    // mc_rtc::log::info("delta {} is\n{}", contactName, delta);

    // bineq vec is torque upper and lower limits (with negative lower limits) + delta (inertia, coriolis, gravity...)
    Eigen::VectorXd bineq(Aineq.rows());

    bineq.segment(0, jacSize) = delta - lowerTorqueLimsNoFloatingBase;
    bineq.segment(jacSize, jacSize) = (-1.0) * (delta - upperTorqueLimsNoFloatingBase);
    // mc_rtc::log::info("b ineq:\n{}", bineq.transpose());

    // Create a half space from every inequality
    for(auto i = 0; i < jacSize * 2; i++)
    {
      // Add half space only if row is not null, ie only if this dof plays into the contact force (reduces nb of planes
      // to simplify in polytope)
      // XXX POLITOPIX: since politopix convention seems to be inverted for Hrep inequalities
      // i.e. they check that a0 + a1*x1 + a2*x2 ... >= 0 to belong inside
      // we push inverted normals
      if(!Aineq.row(i).isZero())
      {
        boost::shared_ptr<HalfSpace_Rn> hs(new HalfSpace_Rn(dim));
        boost::numeric::ublas::vector<double> coefficients(3);
        // Setting coefficients as force elements of the ineq matrix (3 last columns)
        coefficients.insert_element(0, Aineq.coeff(i, 3));
        coefficients.insert_element(1, Aineq.coeff(i, 4));
        coefficients.insert_element(2, Aineq.coeff(i, 5));
        hs->setCoefficients(-coefficients);
        hs->setConstant(bineq.coeff(i));
        newPoly->addHalfSpace(hs);
      }
    }
    auto start_DD = mc_rtc::clock::now();
    // Compute double description from half spaces (not generators -> truncation with bounding box)
    auto result = politopixAPI::computeDoubleDescriptionWithoutCheck(newPoly, 5000);
    // mc_rtc::log::info("DD force poly for {} finished with {} gens and {} hs", contactName,
    //                   newPoly->numberOfGenerators(), newPoly->numberOfHalfSpaces());
    mc_rtc::duration_ms end_DD = mc_rtc::clock::now() - start_DD;
    // mc_rtc::log::info("time to run force poly DD : {}ms", end_DD.count());
  }

  /***************************************************/
  // Vrep version

  else
  {
    // Operational space inertia matrix (OSIM) is (J*M^-1*J^T)^-1
    Eigen::MatrixXd OSIM = (contactFullJac * inertiaMat.inverse() * contactFullJacT).inverse();
    // contactFullJac dim is (40, 6)
    // inertiaMat dim is (40, 40)
    // contactFullJacT dim is (6, 40)
    // OSIM dim is (6, 6)

    // Dynamically consistant inverse jacobian (DCIJ) is M^-1*J^T*OSIM
    // Transposed is equal to OSIM*J*M^-1 because OSIM and joint space inertia matrix are symmetric
    Eigen::MatrixXd DCIJ_T = OSIM * contactFullJac * inertiaMat.inverse();
    // DCIJ_T dims are (6, 40)

    // mc_rtc::log::info("contactFullJac transposed is :\n{}", contactFullJac.transpose());
    /* vertices force polytope are gotten with torque limits vertices
    f = -(JM^-1J^T)^-1JM^-1 * tau + (JM^-1J^T)^-1JM^-1 * C - (JM^-1J^T)^-1*\dot(J)*\dot(q)
    --> f = -DCIJ_T * tau + DCIJ_T * C - OSIM * \dot(J) * \dot(q)
    */

    Eigen::VectorXd constant = DCIJ_T * coriolisVec - OSIM * contactFullJacDot * qdot;

    // Now loop between every combination of min and max torques to get the corresponding 6D wrench
    // Then take only last 3 elements to get the force vertex and add it as generator

    // For the combinations: use only the selected torque extrema (from the actuators of the limb) (2^n combinations and
    // n = nb of actuators)
    auto selectedUpperTorqueLims = actuatorSelectionVector.cwiseProduct(upperTorqueLims);
    auto selectedLowerTorqueLims = actuatorSelectionVector.cwiseProduct(lowerTorqueLims);
    // mc_rtc::log::info("upper torque lims {} are \n{}", contactName, selectedUpperTorqueLims);
    // mc_rtc::log::info("lower torque lims {} are \n{}", contactName, selectedLowerTorqueLims);

    // Lambda to generate combinations
    auto generateCombinations = [](const Eigen::VectorXd & minVec,
                                   const Eigen::VectorXd & maxVec) -> std::vector<Eigen::VectorXd>
    {
      int n = minVec.size();

      // Combine only selected elements : others are given minVec values (I assume they are 0)
      std::vector<int> variableIndices;
      Eigen::VectorXd fixedValues(n);
      for(int i = 0; i < n; i++)
      {
        if(minVec[i] == maxVec[i])
        {
          fixedValues[i] = minVec[i];
        }
        else
        {
          // These are the elements that need to be combined
          variableIndices.push_back(i);
        }
      }

      int nbActuatorsPlaying = variableIndices.size();
      int nbOfCombinations = pow(2, nbActuatorsPlaying);
      std::vector<Eigen::VectorXd> combinations;

      for(int i = 0; i < nbOfCombinations; i++)
      {
        // Start new combination by putting fixed values (variable values are still not set)
        Eigen::VectorXd combination = fixedValues;

        for(int j = 0; j < nbActuatorsPlaying; j++)
        {
          int index = variableIndices[j];
          // Make a combination using either minVec[index] or maxVec[index] depending on j-th bit of i
          combination[index] = (i & (1 << j)) ? maxVec[index] : minVec[index];
        }
        combinations.push_back(combination);
      }
      return combinations;
    };

    // Generate all torque extrema combinations
    std::vector<Eigen::VectorXd> torqueExtremaCombinations =
        generateCombinations(selectedLowerTorqueLims, selectedUpperTorqueLims);

    // mc_rtc::log::info("There are {} combinations", torqueExtremaCombinations.size());

    // Compute resulting force for every combination and add it as vertex for polytope before DD
    int itCounter = 0;

    auto start_Qhull = mc_rtc::clock::now();
    std::vector<double> qhullVect;
    for(const auto & torqueCombination : torqueExtremaCombinations)
    {
      // mc_rtc::log::info("torque combination vector for {} is \n{}", contactName, torqueCombination);
      Eigen::VectorXd wrenchVertex = -DCIJ_T * torqueCombination + constant;
      // mc_rtc::log::info("Wrench vertex {} for {} is {}", itCounter, contactName,
      //                   wrenchVertex.segment(3, 3).transpose());

      // Insert force elements (3, 4, 5 of wrench 6D vector)

      qhullVect.emplace_back(wrenchVertex.coeff(3));
      qhullVect.emplace_back(wrenchVertex.coeff(4));
      qhullVect.emplace_back(wrenchVertex.coeff(5));

      itCounter++;
    }

    computeQhullHrep(qhullVect, newPoly);
    mc_rtc::duration_ms end_Qhull = mc_rtc::clock::now() - start_Qhull;
    // mc_rtc::log::info("time to run force poly Qhull {} : {}ms", contactName, end_Qhull.count());

    auto start_DD = mc_rtc::clock::now();

    // Now newPoly is a polytope with only the H-rep, run DD on it (for intersection)
    // XXX Now keeping Hrep only, see buildFeasiblePolytopeFromContact for logic

    // politopixAPI::computeDoubleDescriptionWithoutCheck(newPoly, 3000);

    // mc_rtc::log::info("{} has {} gens, {} hs after DD", contactName, newPoly->numberOfGenerators(),
    // newPoly->numberOfHalfSpaces());

    // for(int vertex = 0; vertex < newPoly->numberOfGenerators(); vertex++)
    // {
    //   Eigen::Vector3d thisVertex(newPoly->getGenerator(vertex)->getCoordinate(0),
    //                              newPoly->getGenerator(vertex)->getCoordinate(1),
    //                              newPoly->getGenerator(vertex)->getCoordinate(2));
    //   mc_rtc::log::info("Vertex {} of {} at {}", vertex, contactName, thisVertex.transpose());
    // }

    // mc_rtc::log::info("Topo ok {}? {}", contactName, newPoly->checkTopologyAndGeometry());

    mc_rtc::duration_ms end_DD = mc_rtc::clock::now() - start_DD;
    // mc_rtc::log::info("time to run force poly DD {} : {}ms", contactName, end_DD.count());
  }

  // Scaling force polytope before intersection with friction cone
  TopGeomTools::scalingFactor(newPoly, forceScalingFactor);

  // mc_rtc::log::info("Force polytope of {} is ok? {}", contactName, checkGravityCenterInPolytope(newPoly));
  // lock poly mutex, then reset poly pointer to newly computed poly
  std::lock_guard<std::mutex> lock(forcePolyMutex);
  actuationPolytope.reset();
  actuationPolytope = newPoly;
}

void DynamicPolytope::computeQhullHrep(std::vector<double> & points, boost::shared_ptr<Polytope_Rn> & polytope)
{
  // Memory management vars
  int curlong, totlong;
  // Initialize local Qhull environment + pointer to it
  qhT qh_qh;
  qhT * qh = &qh_qh;
  if(!qh)
  {
    mc_rtc::log::error("Qhull malloc error.");
    return;
  }
  QHULL_LIB_CHECK
  // QHull init to be thread safe
  qh_zero(qh, nullptr);

  int dim = 3;
  int numPoints = points.size() / dim;

  // Compute convex envelope (string arg for cpp binding)
  if(qh_new_qhull(qh, dim, numPoints, points.data(), 0, "qhull", nullptr, nullptr))
  {
    mc_rtc::log::error("Qhull: convex envelope error.");
    qh_freeqhull(qh, !qh_ALL);
    qh_memfreeshort(qh, &curlong, &totlong);
    return;
  }

  facetT * facet;
  boost::shared_ptr<HalfSpace_Rn> HS;
  // This macro iterates the qhull environment facet list on the declared facetT
  FORALLfacets
  {
    // Qhull convention is a0 + a1*x1 + a2*x2+ ... = 0
    // This means offset should be negative for inequality
    if(!facet->normal) continue;
    HS.reset(new HalfSpace_Rn(3));
    // Getting all 3 coefficients of the facet normal
    // (negative for politopix convention)
    for(int coord_count = 0; coord_count < 3; coord_count++)
    {
      HS->setCoefficient(coord_count, -facet->normal[coord_count]);
    }
    // Setting offset (negative from qhull convention)
    HS->setConstant(-facet->offset);
    polytope->addHalfSpace(HS);
  }

  // Memory free:
  // Long memory
  qh_freeqhull(qh, !qh_ALL);
  // Short memory and memory allocator
  qh_memfreeshort(qh, &curlong, &totlong);
  if(curlong || totlong)
  {
    mc_rtc::log::error("Qhull memory: did not free {} bytes of long memory ({} pieces)", totlong, curlong);
  }
}

void DynamicPolytope::buildFeasiblePolytopeFromContact(const std::string contactName,
                                                       const mc_rbdyn::Robot & robot,
                                                       const sva::PTransformd X_r1_r2,
                                                       int numberOfFrictionSides,
                                                       double forceScalingFactor,
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
  buildActuationPolytopeFromContact(contactName, robot, actuationPolytope, forcePolyMutex, forceScalingFactor);
  timers.dt_forcePolytope = mc_rtc::clock::now() - start_forcePoly;

  auto start_frictionCone = mc_rtc::clock::now();
  buildFrictionConeFromContactWithHrep(numberOfFrictionSides, X_r1_r2, frictionCone, frictionConeMutex, frictionCoeff);
  timers.dt_frictionCone = mc_rtc::clock::now() - start_frictionCone;

  auto start_intersection = mc_rtc::clock::now();
  // Compute intersection
  std::lock_guard<mutex> lockFriction(frictionConeMutex);
  std::lock_guard<mutex> lockForce(forcePolyMutex);

  // intersect friction cone planes with force polytope into friction cone object
  // Two possible ways: either we don't run the DD on the force poly at first, just compute force planes, compute
  // friction planes, combine them and then run a DD to get intersection: less redundant operations but the DD runs on a
  // big polytope
  // Second way is to run the DD on the force poly as expected then intersect with friction planes. Redundant operations
  // because kind of DD twice but intersection runs only on the newly added planes so in the end almost the same.

  // mc_rtc::log::info("Friction cone {} before intersection: {} hs and {} gens", contactName,
  //                   frictionCones_.at(contactName)->numberOfHalfSpaces(),
  //                   frictionCones_.at(contactName)->numberOfGenerators());
  // mc_rtc::log::info("Force poly {} before intersection: {} hs and {} gens", contactName,
  //                   forcePolytopes_.at(contactName)->numberOfHalfSpaces(),
  //                   forcePolytopes_.at(contactName)->numberOfGenerators());

  // Adding friction cone planes to force planes before running DD

  constIteratorOfListOfGeometricObjects<boost::shared_ptr<HalfSpace_Rn>> iteHSB(frictionCone->getListOfHalfSpaces());
  for(iteHSB.begin(); iteHSB.end() != true; iteHSB.next())
  {
    actuationPolytope->addHalfSpace(iteHSB.current());
  }
  politopixAPI::computeDoubleDescriptionWithoutCheck(actuationPolytope, 3000);

  // politopixAPI::computeIntersectionWithoutCheck(actuationPolytope, frictionCone);

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
    // dummy value: should need difference between controlled surface and the target in the contact pair
    sva::PTransformd X_r1_r2 = sva::PTransformd::Identity();
    if(!withMoments_)
    {
      buildFrictionConeFromContactWithHrep(nbFrictionSides, X_r1_r2, frictionCones_.at(contactName),
                                           getContactMutex(frictionConesMutexes_, contactName), frictionCoeff);
    }
    else
    {
      // find limits of contact area for moment limits
      double newContactHalfLength;
      double newContactHalfWidth;
      findHalfWidthLength(robot.surface(contactName), newContactHalfWidth, newContactHalfLength);
      std::pair<std::pair<double, double>, sva::PTransformd> newContact(
          std::pair<double, double>(newContactHalfLength, newContactHalfWidth), robot.surfacePose(contactName));

      // TODO thread moments versions as well: add moment mutex + put mutexes as arguments
      buildWrenchConeFromContact(nbFrictionSides, newContact, frictionCones_.at(contactName),
                                 frictionConesMoments_.at(contactName), frictionCoeff, maxForce, CoM);
    }
  }
}

void DynamicPolytope::computeForcePolyFromContactSet(const mc_rbdyn::Robot & robot)
{
  double forceScalingFactor = 1;
  for(const auto contactName : activeContacts_)
  {
    buildActuationPolytopeFromContact(contactName, robot, forcePolytopes_.at(contactName),
                                      getContactMutex(forcePolyMutexes_, contactName), forceScalingFactor);
  }
}

void DynamicPolytope::computeFeasibleForcesFromContactSet(const mc_rbdyn::Robot & robot)
{
  Rn::setDimension(3);
  std::lock_guard<std::mutex> lock(contactSetMutex_);
  for(const auto & contactName : activeContacts_)
  {
    // launching contact computation
    std::lock_guard<std::mutex> lockThreadMap(feasiblePolytopesThreadsMutex_);
    feasiblePolytopesThreads_.emplace(
        contactName,
        std::thread(
            &DynamicPolytope::buildFeasiblePolytopeFromContact, this, contactName, std::ref(robot),
            refContactTransforms_.at(contactName), std::ref(numbersOfFrictionSides_.at(contactName)),
            std::ref(forceScalingFactors_.at(contactName)), std::ref(frictionCoefficients_.at(contactName)),
            std::ref(frictionCones_.at(contactName)), std::ref(getContactMutex(frictionConesMutexes_, contactName)),
            std::ref(forcePolytopes_.at(contactName)), std::ref(getContactMutex(forcePolyMutexes_, contactName)),
            std::ref(contactsTimers_.at(contactName))));
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
  // mc_rtc::log::info("Starting to add polytopes for mink sum");
  for(const auto & active : activeContacts_)
  {
    boost::shared_ptr<Polytope_Rn> newContactPoly(new Polytope_Rn());
    // The friction + force intersection was overwritten in the force polytopes
    // Locking polytope mutex
    // Copying the contact polytopes to use for mink computation
    getContactMutex(forcePolyMutexes_, active).lock();
    politopixAPI::copyPolytope(forcePolytopes_.at(active), newContactPoly);
    // mc_rtc::log::info("Adding poly with {} gens and {} hs", newContactPoly->numberOfGenerators(),
    //                   newContactPoly->numberOfHalfSpaces());
    getContactMutex(forcePolyMutexes_, active).unlock();

    // After copying the contact frame polytope, rotate it to world frame before minkowski sum
    // Vectors are force in contact frame so X_contact_f
    // We want X_0_f = X_contact_f * X_0_contact
    auto X_0_contact = robot_.surfacePose(active);
    X_0_contact.translation() = Eigen::Vector3d::Zero();

    // Rotating generators
    constIteratorOfListOfGeometricObjects<boost::shared_ptr<Generator_Rn>> iterGen(
        newContactPoly->getListOfGenerators());
    for(iterGen.begin(); iterGen.end() != true; iterGen.next())
    {
      Eigen::Vector3d vect(iterGen.current()->getCoordinate(0), iterGen.current()->getCoordinate(1),
                           iterGen.current()->getCoordinate(2));
      vect = (sva::PTransformd(vect) * X_0_contact).translation();
      iterGen.current()->setCoordinate(0, vect.x());
      iterGen.current()->setCoordinate(1, vect.y());
      iterGen.current()->setCoordinate(2, vect.z());
    }
    // Rotating halfspaces normals
    constIteratorOfListOfGeometricObjects<boost::shared_ptr<HalfSpace_Rn>> iterHS(
        newContactPoly->getListOfHalfSpaces());
    for(iterHS.begin(); iterHS.end() != true; iterHS.next())
    {
      Eigen::Vector3d normal(iterHS.current()->getCoefficient(0), iterHS.current()->getCoefficient(1),
                             iterHS.current()->getCoefficient(2));
      normal = (sva::PTransformd(normal) * X_0_contact).translation();
      iterHS.current()->setCoefficient(0, normal.x());
      iterHS.current()->setCoefficient(1, normal.y());
      iterHS.current()->setCoefficient(2, normal.z());
    }

    polytopesForces.emplace_back(newContactPoly);
    if(withMoments_)
    {
      polytopesMoments.emplace_back(frictionConesMoments_.at(active));
    }
  }
  contactSetMutex_.unlock();

  try
  {
    if(!polytopesForces.empty())
    {
      MinkowskiSum Mink(polytopesForces, newForcePoly);
    }
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
  std::lock_guard<std::mutex> lockCWC(CWCMutex_);
  TopGeomTools::translate(CWCForces_, deltaZVector);
  std::lock_guard<std::mutex> lockZeroMoment(zeroMomentMutex_);
  TopGeomTools::translate(zeroMomentRegion_, deltaZVector);
}

void DynamicPolytope::computeZMPRegion(Eigen::Vector3d comPosition)
{
  // XXX dummy zone for now: convex area formed by the polygon envelope of feet + com position
  int dim = 3;
  boost::shared_ptr<Polytope_Rn> newZMPPoly(new Polytope_Rn());

  std::vector<double> qhullVect;
  const std::string contactTest = "LeftFoot";
  for(const auto & contact : activeContacts_)
  {
    // Important ! points for surfaces give the points coordinates from the parent link, not from the
    // surface origin, so X_b_p and not X_s_p
    // This means to get world frame surface points X_0_p = X_s_p * X_0_s we need X_s_p
    // X_s_p = X_b_p * X_s_b = X_b_p * X_b_s.inv()
    auto cPoints = robot_.surface(contact).points();
    for(auto cPoint : cPoints)
    {
      cPoint = cPoint * robot_.surface(contact).X_b_s().inv() * robot_.surface(contact).X_0_s(robot_);
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

  computeQhullHrep(qhullVect, newZMPPoly);

  politopixAPI::computeDoubleDescriptionWithoutCheck(newZMPPoly, 3000);

  checkAllHSInternal("ZMP", newZMPPoly);

  std::lock_guard<std::mutex> lock(ZMPMutex_);
  zmpRegion_ = newZMPPoly;
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
  if(zmpThread_.joinable())
  {
    zmpThread_.join();
  }
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

void DynamicPolytope::checkAllHSInternal(const std::string & polyName, boost::shared_ptr<Polytope_Rn> & polytope)
{
  // Checking if there are internal halfspaces
  constIteratorOfListOfGeometricObjects<boost::shared_ptr<HalfSpace_Rn>> HSiter(polytope->getListOfHalfSpaces());
  for(HSiter.begin(); HSiter.end() != true; HSiter.next())
  {
    // Checking if every vertex is inside of this HS
    constIteratorOfListOfGeometricObjects<boost::shared_ptr<Generator_Rn>> GNiter(polytope->getListOfGenerators());
    for(GNiter.begin(); GNiter.end() != true; GNiter.next())
    {
      double scalarProduct =
          std::inner_product(GNiter.current()->begin(), GNiter.current()->end(), HSiter.current()->begin(), 0.);
      if(scalarProduct + HSiter.current()->getConstant() < -1.e-04)
      {
        mc_rtc::log::critical("{} polytope has an internal HS (error {})", polyName,
                              scalarProduct + HSiter.current()->getConstant());
      }
    }
  }
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
    auto contactPose = robot_.surfacePose(contact);
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
  int nbOfPlanes = polytope->numberOfHalfSpaces();
  Normals.resize(nbOfPlanes, 3);
  Offsets.resize(nbOfPlanes);

  for(int halfSpaceIndex = 0; halfSpaceIndex < nbOfPlanes; halfSpaceIndex++)
  {
    const auto halfSpace = polytope->getHalfSpace(halfSpaceIndex);
    // XXX POLITOPIX: since politopix convention seems to be inverted for Hrep inequalities
    // i.e. they check that a0 + a1*x1 + a2*x2 ... >= 0 to belong inside
    // we push inverted normals
    Normals.row(halfSpaceIndex) << -halfSpace->getCoefficient(0), -halfSpace->getCoefficient(1),
        -halfSpace->getCoefficient(2);
    Offsets(halfSpaceIndex) = halfSpace->getConstant();
  }
}

void DynamicPolytope::updateRBDynPolytopes(const Eigen::MatrixX3d & Normals,
                                           const Eigen::VectorXd & Offsets,
                                           mc_rbdyn::Contact & contactRBDyn)
{
  mc_rbdyn::FeasiblePolytope polytope({Normals, Offsets});
  // Update the contact polytope of the desired robot
  if(robot_.robotIndex() == contactRBDyn.r1Index())
  {
    contactRBDyn.feasiblePolytopeR1(polytope);
  }
  else
  {
    contactRBDyn.feasiblePolytopeR2(polytope);
  }
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

void DynamicPolytope::addToGUI(mc_rtc::gui::StateBuilder & gui, double guiScale, std::vector<std::string> category)
{
  config_("gui")("polyhedronForce")("triangle_color", polyForceConfig_.triangle_color);
  config_("gui")("polyhedronForce")("show_triangle", polyForceConfig_.show_triangle);
  config_("gui")("polyhedronForce")("use_triangle_color", polyForceConfig_.use_triangle_color);
  config_("gui")("polyhedronForce")("edges", polyForceConfig_.edge_config);
  config_("gui")("polyhedronForce")("show_edges", polyForceConfig_.show_edges);
  config_("gui")("polyhedronForce")("fixed_edge_color", polyForceConfig_.fixed_edge_color);
  config_("gui")("polyhedronForce")("vertices")("color", polyForceConfig_.vertices_config.color);
  config_("gui")("polyhedronForce")("vertices")("scale", polyForceConfig_.vertices_config.scale);
  config_("gui")("polyhedronForce")("show_vertices", polyForceConfig_.show_vertices);
  config_("gui")("polyhedronForce")("fixed_vertices_color", polyForceConfig_.fixed_vertices_color);

  config_("gui")("polyhedronMoment")("triangle_color", polyMomentConfig_.triangle_color);
  config_("gui")("polyhedronMoment")("show_triangle", polyMomentConfig_.show_triangle);
  config_("gui")("polyhedronMoment")("use_triangle_color", polyMomentConfig_.use_triangle_color);
  config_("gui")("polyhedronMoment")("edges", polyMomentConfig_.edge_config);
  config_("gui")("polyhedronMoment")("show_edges", polyMomentConfig_.show_edges);
  config_("gui")("polyhedronMoment")("fixed_edge_color", polyMomentConfig_.fixed_edge_color);
  config_("gui")("polyhedronMoment")("vertices")("color", polyMomentConfig_.vertices_config.color);
  config_("gui")("polyhedronMoment")("vertices")("scale", polyMomentConfig_.vertices_config.scale);
  config_("gui")("polyhedronMoment")("show_vertices", polyMomentConfig_.show_vertices);
  config_("gui")("polyhedronMoment")("fixed_vertices_color", polyMomentConfig_.fixed_vertices_color);

  config_("gui")("polyhedronZMP")("triangle_color", polyZMPConfig_.triangle_color);
  config_("gui")("polyhedronZMP")("show_triangle", polyZMPConfig_.show_triangle);
  config_("gui")("polyhedronZMP")("use_triangle_color", polyZMPConfig_.use_triangle_color);
  config_("gui")("polyhedronZMP")("edges", polyZMPConfig_.edge_config);
  config_("gui")("polyhedronZMP")("show_edges", polyZMPConfig_.show_edges);
  config_("gui")("polyhedronZMP")("fixed_edge_color", polyZMPConfig_.fixed_edge_color);
  config_("gui")("polyhedronZMP")("vertices")("color", polyZMPConfig_.vertices_config.color);
  config_("gui")("polyhedronZMP")("vertices")("scale", polyZMPConfig_.vertices_config.scale);
  config_("gui")("polyhedronZMP")("show_vertices", polyZMPConfig_.show_vertices);
  config_("gui")("polyhedronZMP")("fixed_vertices_color", polyZMPConfig_.fixed_vertices_color);

  config_("gui")("polyhedronZeroMomentArea")("triangle_color", polyZeroMomentAreaConfig_.triangle_color);
  config_("gui")("polyhedronZeroMomentArea")("show_triangle", polyZeroMomentAreaConfig_.show_triangle);
  config_("gui")("polyhedronZeroMomentArea")("use_triangle_color", polyZeroMomentAreaConfig_.use_triangle_color);
  config_("gui")("polyhedronZeroMomentArea")("edges", polyZeroMomentAreaConfig_.edge_config);
  config_("gui")("polyhedronZeroMomentArea")("show_edges", polyZeroMomentAreaConfig_.show_edges);
  config_("gui")("polyhedronZeroMomentArea")("fixed_edge_color", polyZeroMomentAreaConfig_.fixed_edge_color);
  config_("gui")("polyhedronZeroMomentArea")("vertices")("color", polyZeroMomentAreaConfig_.vertices_config.color);
  config_("gui")("polyhedronZeroMomentArea")("vertices")("scale", polyZeroMomentAreaConfig_.vertices_config.scale);
  config_("gui")("polyhedronZeroMomentArea")("show_vertices", polyZeroMomentAreaConfig_.show_vertices);
  config_("gui")("polyhedronZeroMomentArea")("fixed_vertices_color", polyZeroMomentAreaConfig_.fixed_vertices_color);

  guiScale_ = guiScale;
  category.push_back(name_);
  auto coeffsCat = category;
  coeffsCat.push_back("Coefficients");
  auto contactsCat = category;
  contactsCat.push_back("Contact Polytopes");
  auto CWCCat = category;
  CWCCat.push_back("Contact Wrench Cone");

  gui.addElement(
      this, category,
      mc_rtc::gui::Checkbox(
          "Compute explicit regions", [this]() { return computeRegions_; },
          [this]() { computeRegions_ = !computeRegions_; }),
      mc_rtc::gui::Checkbox(
          "Compute force poly from Hrep", [this]() { return HrepMode_; }, [this]() { HrepMode_ = !HrepMode_; }));

  for(const auto contact : possibleContacts_)
  {
    gui.addElement(this, coeffsCat,
                   mc_rtc::gui::NumberSlider(
                       fmt::format(contact + " force alpha [0.001-1]"),
                       [this, contact]() { return getForceScalingFactor(contact); },
                       // scale min at 0.001 instead of 0 to avoid handling zero polytope
                       [this, contact](double scale) { getForceScalingFactor(contact) = scale; }, 0.001, 1.0),
                   mc_rtc::gui::NumberInput(
                       fmt::format(contact + " friction coefficient"),
                       [this, contact]() { return getFrictionCoeff(contact); },
                       [this, contact](double frictionCoeff) { getFrictionCoeff(contact) = frictionCoeff; }),
                   mc_rtc::gui::IntegerInput(
                       fmt::format(contact + " number of friction sides"),
                       [this, contact]() { return getNbOfFrictionSides(contact); },
                       [this, contact](int nbFrictionSides) { getNbOfFrictionSides(contact) = nbFrictionSides; }));

    gui.addElement(this, contactsCat,
                   mc_rtc::gui::Polyhedron(fmt::format(contact + " frictions"), polyForceConfig_,
                                           [this, contact]() { return getFrictionConesTriangles(contact); }),
                   mc_rtc::gui::Polyhedron(fmt::format(contact + " forces"), polyMomentConfig_,
                                           [this, contact]() { return getForcePolyTriangles(contact); })/*,
                   mc_rtc::gui::Polyhedron(fmt::format(contact + " moments"), polyMomentConfig_,
                                           [this, contact]() { return getContactMomentTriangles(contact); })*/);
  }

  gui.addElement(
      this, CWCCat,
      mc_rtc::gui::Polyhedron("CWC forces", polyForceConfig_, [this]() { return getCWCForceTriangles(); }),
      /*mc_rtc::gui::Polyhedron("CWC moments", polyMomentConfig_, [this]() { return getCWCMomentTriangles(); }),*/
      mc_rtc::gui::Polyhedron("ZMP area", polyZMPConfig_, [this]() { return getZMPTriangles(); }),
      mc_rtc::gui::Polyhedron("Zero moment region", polyZeroMomentAreaConfig_,
                              [this]() { return getZeroMomentTriangles(); }));
}
