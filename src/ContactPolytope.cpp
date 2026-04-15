#include <RBDyn/Coriolis.h>
#include <mc_dynamic_polytopes/ContactPolytope.h>

namespace mc_dynamic_polytopes
{

void ContactPolytopeJob::buildActuationPolytopeFromContact(const ContactPolytopeInput & input,
                                                           boost::shared_ptr<Polytope_Rn> & actuationPolytope,
                                                           double forceScalingFactor,
                                                           unsigned int dim)
{
  boost::shared_ptr<Polytope_Rn> newPoly(new Polytope_Rn());

  const int n_var = 6;
  // Removing underactuated dofs
  const int jacSize = input.mb.nrDof() - 6;

  // TODO: pass this option in inputs
  // This decides if we compute the polytope using only the force jacobian (true) or using the 6D jacobian and taking
  // only the associated force part for visualization (false)
  // Setting it to true should be incompatible with the withMoments option (that sets dim to 6)
  bool forceOnly = false;

  // Computing dynamics terms of the equation of motion
  rbd::ForwardDynamics forwardDyn(input.mb);
  // Inertia
  forwardDyn.computeH(input.mb, input.mbc);
  Eigen::MatrixXd inertiaMat = forwardDyn.H();

  // Coriolis
  // XXX should the coriolis matrix be used or the C vector of the forward dynamics?
  rbd::Coriolis coriolis(input.mb);
  Eigen::MatrixXd coriolisMat = coriolis.coriolis(input.mb, input.mbc);

  // XXX this version contains both gravity and external forces, we use this one directly
  forwardDyn.computeC(input.mb, input.mbc);
  Eigen::VectorXd coriolisVec = forwardDyn.C();

  // Correct jacobian building needs names because we need to link the surface to a body and compute the jac of the body
  // + transform between surface and the body
  const mc_rbdyn::Surface & surface = *input.surface;
  std::string bodyName = surface.bodyName();
  const auto X_b_s = surface.X_b_s();

  // Building jacobian of the frame point on body in world frame
  rbd::Jacobian jac(input.mb, bodyName, X_b_s.translation());

  // Dense jacobian in body frame
  Eigen::MatrixXd denseJac = jac.bodyJacobian(input.mb, input.mbc);
  // Dense jac dot in body frame
  Eigen::MatrixXd denseJacDot = jac.bodyJacobianDot(input.mb, input.mbc);
  // mc_rtc::log::info("Dense Jacobian: \n{}", denseJac);
  // mc_rtc::log::info("Dense Jacobian {} is of dimensions {}", contactName, denseJac.rows(), denseJac.cols());

  // Allocate then fill sparse jacobian
  Eigen::MatrixXd globalFullJac = Eigen::MatrixXd::Zero(6, input.mb.nrDof());
  jac.fullJacobian(input.mb, denseJac, globalFullJac);
  Eigen::MatrixXd globalFullJacDot = Eigen::MatrixXd::Zero(6, input.mb.nrDof());
  jac.fullJacobian(input.mb, denseJacDot, globalFullJacDot);
  // mc_rtc::log::info("Sparse Jacobian: \n{}", globalFullJac);

  // rotate from body to contact
  Eigen::MatrixXd contactFullJac6D = X_b_s.matrix() * globalFullJac;
  // XXX to check ! simple rotation if jacDot matrix?
  Eigen::MatrixXd contactFullJacDot6D = X_b_s.matrix() * globalFullJac;

  // If we want force only, just take the bottom 3 rows of the jacs, otherwise full
  Eigen::MatrixXd contactFullJac;
  Eigen::MatrixXd contactFullJacDot;
  if(forceOnly)
  {
    contactFullJac = contactFullJac6D.bottomRows(3);
    contactFullJacDot = contactFullJacDot6D.bottomRows(3);
  }
  else
  {
    contactFullJac = contactFullJac6D;
    contactFullJacDot = contactFullJacDot6D;
  }

  // transpose for correct calculation
  const Eigen::MatrixXd contactFullJacT = contactFullJac.transpose();
  // mc_rtc::log::info("Sparse jacobian transposed and rotated to surface frame: \n{}", contactFullJacT);

  // Torque limits vectors
  Eigen::VectorXd upperTorqueLims = rbd::dofToVector(input.mb, input.tu);
  Eigen::VectorXd lowerTorqueLims = rbd::dofToVector(input.mb, input.tl);
  Eigen::VectorXd upperTorqueLimsNoFloatingBase = rbd::dofToVector(input.mb, input.tu).segment(6, jacSize);
  Eigen::VectorXd lowerTorqueLimsNoFloatingBase = rbd::dofToVector(input.mb, input.tl).segment(6, jacSize);

  // Joint vectors
  const Eigen::VectorXd qdot = rbd::dofToVector(input.mb, input.mbc.alpha);
  const Eigen::VectorXd qddot = rbd::dofToVector(input.mb, input.mbc.alphaD);

  // Selection vector for contact by using the contact jacobian: select this dof if jacobian of contact is not zero
  Eigen::VectorXd actuatorSelectionVector = Eigen::VectorXd::Zero(input.mb.nrDof());
  int numberOfActuatorsPlaying = 0;

  for(int i = 0; i < input.mb.nrDof(); i++)
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

    Eigen::VectorXd delta = inertiaMat.bottomRows(jacSize).leftCols(6) * input.accW.vector()
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
        boost::numeric::ublas::vector<double> coefficients(dim);
        // Setting coefficients as force or wrenches elements of the ineq matrix (3 last columns or all 6)
        if(dim == 3)
        {
          coefficients.insert_element(0, Aineq.coeff(i, 3));
          coefficients.insert_element(1, Aineq.coeff(i, 4));
          coefficients.insert_element(2, Aineq.coeff(i, 5));
        }
        if(dim == 6)
        {
          coefficients.insert_element(0, Aineq.coeff(i, 0));
          coefficients.insert_element(1, Aineq.coeff(i, 1));
          coefficients.insert_element(2, Aineq.coeff(i, 2));
          coefficients.insert_element(3, Aineq.coeff(i, 3));
          coefficients.insert_element(4, Aineq.coeff(i, 4));
          coefficients.insert_element(5, Aineq.coeff(i, 5));
        }

        hs->setCoefficients(-coefficients);
        hs->setConstant(bineq.coeff(i));
        newPoly->addHalfSpace(hs);
      }
    }
    auto start_DD = mc_rtc::clock::now();
    // Compute double description from half spaces (not generators -> truncation with bounding box)
    politopixAPI::computeDoubleDescriptionWithoutCheck(newPoly, 5000);
    // mc_rtc::log::info("DD force poly for {} finished with {} gens and {} hs", contactName,
    //                   newPoly->numberOfGenerators(), newPoly->numberOfHalfSpaces());
    contactTimers.dt_double_description = mc_rtc::elapsed_ms_count(start_DD);
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
    f = -(JM^-1J^T)^-1JM^-1 * tau + (JM^-1J^T)^-1JM^-1 * C - (JM^-1J^T)^-1*\dot(J)*\dot(q) + (JM^-1J^T)^-1 * \ddot(p)
    --> f = -DCIJ_T * tau + DCIJ_T * C - OSIM * \dot(J) * \dot(q) + OSIM \ddot(p)
    with p being the desired task space acceleration of the contact
    */

    // TODO: For now we choose zero, find a way to pass it to job, need one per contact
    Eigen::VectorXd contactAcc;
    if(forceOnly)
    {
      contactAcc = Eigen::Vector3d::Zero();
    }
    else
    {
      contactAcc = Eigen::Vector6d::Zero();
    }

    Eigen::VectorXd constant = DCIJ_T * coriolisVec - OSIM * contactFullJacDot * qdot + OSIM * contactAcc;

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
      // mc_rtc::log::info("Wrench vertex {} for {} is {}", itCounter, input.contactName,
      //                   wrenchVertex.transpose());

      // Insert force elements (3, 4, 5 of wrench 6D vector) or all wrench
      if(dim == 6)
      {
        qhullVect.emplace_back(wrenchVertex.coeff(0));
        qhullVect.emplace_back(wrenchVertex.coeff(1));
        qhullVect.emplace_back(wrenchVertex.coeff(2));
        qhullVect.emplace_back(wrenchVertex.coeff(3));
        qhullVect.emplace_back(wrenchVertex.coeff(4));
        qhullVect.emplace_back(wrenchVertex.coeff(5));
      }
      else // dim is 3, check whether we compute the force jac only or the force part of the full jac
      {
        if(forceOnly)
        {
          qhullVect.emplace_back(wrenchVertex.coeff(0));
          qhullVect.emplace_back(wrenchVertex.coeff(1));
          qhullVect.emplace_back(wrenchVertex.coeff(2));
        }
        else
        {
          qhullVect.emplace_back(wrenchVertex.coeff(3));
          qhullVect.emplace_back(wrenchVertex.coeff(4));
          qhullVect.emplace_back(wrenchVertex.coeff(5));
        }
      }

      itCounter++;
    }

    computeQhullHrep(qhullVect, newPoly, dim);
    contactTimers.dt_qhull = mc_rtc::elapsed_ms_count(start_Qhull);
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

    mc_rtc::duration_ms end_DD = mc_rtc::elapsed_ms(start_DD);
    // mc_rtc::log::info("time to run force poly DD {} : {}ms", contactName, end_DD.count());
  }

  // Scaling force polytope before intersection with friction cone
  TopGeomTools::scalingFactor(newPoly, forceScalingFactor);

  // mc_rtc::log::info("Force polytope of {} is ok? {}", contactName, checkGravityCenterInPolytope(newPoly));
  // lock poly mutex, then reset poly pointer to newly computed poly

  // XXX: do we need to reset the pointer?
  actuationPolytope.reset();
  actuationPolytope = newPoly;
}

ContactPolytopeResult ContactPolytopeJob::computeJob()
{
  ContactPolytopeResult result;
  auto & frictionCone = result.frictionCone;
  auto & actuationPolytope = result.actuationPolytope;
  const auto X_r1_r2 = input_.refContactTransform; // XXX: double check
  const auto frictionCoeff = input_.frictionCoefficient;

  // sva::PTransformd contactPose = robot.surfacePose(contactName);
  unsigned int dim = Rn::getDimension() == 3 ? 3 : 6;
  auto start_forcePoly = mc_rtc::clock::now();

  // update the correct force polytope in the map
  buildActuationPolytopeFromContact(input_, actuationPolytope, input_.forceScalingFactor, dim);
  contactTimers.dt_forcePolytope = mc_rtc::elapsed_ms_count(start_forcePoly);

  auto start_frictionCone = mc_rtc::clock::now();
  buildFrictionConeFromContactWithHrep(input_.numberOfFrictionSides, X_r1_r2, frictionCone, frictionCoeff, dim);
  contactTimers.dt_frictionCone = mc_rtc::elapsed_ms_count(start_frictionCone);

  auto start_intersection = mc_rtc::clock::now();
  // Compute intersection

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

  // Run DD on friction cones for display (makes it wrong because polyhedral cones should not be bounded)
  if(DDfrictionCones_)
  {
    politopixAPI::computeDoubleDescriptionWithoutCheck(frictionCone, 1000);
  }

  // Adding friction cone planes to force planes before running DD

  if(combineWithFriction_)
  {
    constIteratorOfListOfGeometricObjects<boost::shared_ptr<HalfSpace_Rn>> iteHSB(frictionCone->getListOfHalfSpaces());
    for(iteHSB.begin(); iteHSB.end() != true; iteHSB.next())
    {
      actuationPolytope->addHalfSpace(iteHSB.current());
    }
  }

  if(dim == 3)
  {
    politopixAPI::computeDoubleDescriptionWithoutCheck(actuationPolytope, 5000);
  }

  // politopixAPI::computeIntersectionWithoutCheck(actuationPolytope, frictionCone);

  // mc_rtc::log::info("Friction cone {} after intersection: {} hs and {} gens", contactName,
  //                   frictionCones_.at(contactName)->numberOfHalfSpaces(),
  //                   frictionCones_.at(contactName)->numberOfGenerators());
  // mc_rtc::log::info("Force poly {} after intersection: {} hs and {} gens", contactName,
  //                   forcePolytopes_.at(contactName)->numberOfHalfSpaces(),
  //                   forcePolytopes_.at(contactName)->numberOfGenerators());
  contactTimers.dt_intersection = mc_rtc::elapsed_ms_count(start_intersection);
  contactTimers.dt_contactTotal = mc_rtc::elapsed_ms_count(start_forcePoly);

  updatePlanesMatrixConstraint(result.frictionCone, result.frictionConesPlanes.first,
                               result.frictionConesPlanes.second);
  updatePlanesMatrixConstraint(result.actuationPolytope, result.forcePolyPlanes.first, result.forcePolyPlanes.second);

  // Check if a gui scale exists in datastore from somewhere else, and if so apply it to the polytope (GUI only)
  if(ctl_)
  {
    if(ctl_->datastore().has("Polytopes::GUIScale::" + input_.contactName))
    {
      double scale = ctl_->datastore().call<double>("Polytopes::GUIScale::" + input_.contactName);
      result.updateTrianglesGUIPolytopix(scale * guiScale_, input_.surfacePose);
    }
    else
    {
      result.updateTrianglesGUIPolytopix(guiScale_, input_.surfacePose);
    }
  }
  else
  {
    result.updateTrianglesGUIPolytopix(guiScale_, input_.surfacePose);
  }

  return result;
}
} // namespace mc_dynamic_polytopes
