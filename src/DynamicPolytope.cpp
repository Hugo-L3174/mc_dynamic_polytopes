#include "DynamicPolytope.h"

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
    std::vector<std::array<Eigen::Vector3d, 3UL>> newForceTrianglesArray;
    forceConesTrianglesMap_.emplace(contact, newForceTrianglesArray);
    momentPolytopesTrianglesMap_.emplace(contact, newForceTrianglesArray);

    // init cones maps
    boost::shared_ptr<Polytope_Rn> newCone(new Polytope_Rn());
    frictionCones_.emplace(contact, newCone);

    boost::shared_ptr<Polytope_Rn> newMomentCone(new Polytope_Rn());
    frictionConesMoments_.emplace(contact, newMomentCone);
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
  logger.addLogEntry(prefix + "dt_totalLoop", this, [this]() { return dt_loop_total().count(); });
  logger.addLogEntry(prefix + "dt_computeContactSet", this, [this]() { return dt_contactSet().count(); });
  logger.addLogEntry(prefix + "dt_minkSum", this, [this]() { return dt_minkSum().count(); });
  logger.addLogEntry(prefix + "dt_updatePlanes", this, [this]() { return dt_updatePlanes().count(); });
  logger.addLogEntry(prefix + "dt_guiTriangles", this, [this]() { return dt_guiTriangles().count(); });
  logger.addLogEntry(prefix + "dt_zeroMomentIntersection", this, [this]() { return dt_zeroMomentInter().count(); });
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
    // Step 1.1: launch individual cones calculations in separate threads as they are independant
    computeConesFromContactSet(robot_);

    // Step 1.2: launch ZMP region calculation at the same time, also independant
    zmpThread_ = std::thread(&DynamicPolytope::computeZMPRegion, this, robot_.com());

    // Step 2: wait for finished individual cones to start mink sum of the cones
    for(auto contactName : activeContacts_)
    {
      // join all friction cones threads, ie wait they finished, then erase them from the map
      // mc_rtc::log::info("Joining {} thread and erasing it", contactName);
      frictionConesThreads_.at(contactName).join();
      frictionConesThreads_.erase(contactName);
    }
    dt_compute_contactSet_ = mc_rtc::clock::now() - start_loop;

    // All friction cones computed, start minkowsky sum computation
    auto start_minkSum = mc_rtc::clock::now();
    computeMinkowskySumPolitopix();

    // Step 3: wait for finished mink sum to convert to eCMP region (no thread needed)
    dt_compute_minkSum_ = mc_rtc::clock::now() - start_minkSum;
    computeECMPRegion(robot_.com(), robot_);
    // Update eCMP planes internal variables to be fetched by controller
    auto start_updatePlanes = mc_rtc::clock::now();
    eCMPPlanesMutex_.lock();
    updatePlanesMatrixConstraint(CWCForces_, eCMPNormals_, eCMPOffsets_);
    eCMPPlanesMutex_.unlock();
    dt_update_planes_ = mc_rtc::clock::now() - start_updatePlanes;

    // Step 4: wait for finished ZMP region to start zero moment intersection with eCMP region
    zmpThread_.join();
    auto start_zeroMomentIntersection = mc_rtc::clock::now();
    computeZeroMomentIntersection();
    dt_zeroMoment_intersection_ = mc_rtc::clock::now() - start_zeroMomentIntersection;

    // Update zero moment region planes to be fetched by controller
    zeroMomentPlanesMutex_.lock();
    updatePlanesMatrixConstraint(zeroMomentRegion_, zeroMomentNormals_, zeroMomentOffsets_);
    zeroMomentPlanesMutex_.unlock();

    // Step 5: gui computations
    auto start_guiTriangles = mc_rtc::clock::now();
    updateTrianglesGUIPolitopix();
    dt_compute_guiTriangles_ = mc_rtc::clock::now() - start_guiTriangles;
    dt_loop_total_ = mc_rtc::clock::now() - start_loop;
  }
}

void DynamicPolytope::buildForceConeFromContact(int numberOfFrictionSides,
                                                sva::PTransformd contactSurface,
                                                boost::shared_ptr<Polytope_Rn> & forceCone,
                                                std::mutex & forceConeMutex,
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
  std::lock_guard<std::mutex> lock(forceConeMutex);
  forceCone.reset();
  forceCone = newCone;
  // mc_rtc::log::info("Created cone of dim {} with {} generators", forceCone->dimension(),
  //                   forceCone->numberOfGenerators());
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

void DynamicPolytope::buildActuationPolytopeFromContact(boost::shared_ptr<Polytope_Rn> & actuationPolytope)
{
  // need to find the branches belonging to the contacts

  actuationPolytope->reset();
}

void DynamicPolytope::computeConesFromContactSet(const mc_rbdyn::Robot & robot)
{
  Rn::setDimension(3);
  auto frictionCoeff = 0.7;
  auto nbFrictionSides = 5;
  auto maxForce = 180.;
  auto CoM = robot.com();

  for(const auto contactName : activeContacts_)
  {
    sva::PTransformd contactPose = robot.surfacePose(contactName);
    if(!withMoments_)
    {
      // update the correct cone in the map
      // launching thread and emplacing it in the threads map
      frictionConesThreads_.emplace(contactName,
                                    std::thread(&DynamicPolytope::buildForceConeFromContact, this, nbFrictionSides,
                                                contactPose, std::ref(frictionCones_.at(contactName)),
                                                std::ref(getContactMutex(frictionConesMutexes_, contactName)),
                                                frictionCoeff, maxForce));
#ifndef WIN32
      // Lower thread priority so that it has a lesser priority than the real time thread
      auto th_handle = frictionConesThreads_.at(contactName).native_handle();
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

void DynamicPolytope::computeMinkowskySumPolitopix()
{
  CWCForces_->reset();
  CWCMoments_->reset();
  // putting it in vector form for library function
  std::vector<boost::shared_ptr<Polytope_Rn>> polytopesForces;
  std::vector<boost::shared_ptr<Polytope_Rn>> polytopesMoments;
  for(auto active : activeContacts_)
  {
    polytopesForces.emplace_back(frictionCones_.at(active));
    if(withMoments_)
    {
      polytopesMoments.emplace_back(frictionConesMoments_.at(active));
    }
  }

  MinkowskiSum Mink(polytopesForces, CWCForces_);
  // mc_rtc::log::info("CWCForces_ has {} generators and {} facets", CWCForces_->numberOfGenerators(),
  // CWCForces_->numberOfHalfSpaces());
  if(withMoments_)
  {
    MinkowskiSum Mink(polytopesMoments, CWCMoments_);
    // mc_rtc::log::info("CWCMoments_ has {} generators and {} facets", CWCMoments_->numberOfGenerators(),
    // CWCMoments_->numberOfHalfSpaces());
  }
}

void DynamicPolytope::computeECMPRegion(Eigen::Vector3d comPosition, const mc_rbdyn::Robot & robot)
{
  // First we scale the force polytope according to the eCMP expression:
  // eCMP = c - sumF/(m*(g/Dz)) --> eCMP = c - sumF*(Dz/mg)
  // scale force polytope by - Deltaz/mg
  double scale = -comPosition.z() / (robot.mass() * 9.81);
  bool ok = TopGeomTools::scalingFactor(CWCForces_, scale);
  // CWCForces_->negate();
  // then translate it from origin to the robot CoM
  boost::numeric::ublas::vector<double> CoM(3);
  CoM[0] = comPosition.x();
  CoM[1] = comPosition.y();
  CoM[2] = comPosition.z();
  TopGeomTools::translate(CWCForces_, CoM);
  // add Delta Z for VRP
  // CoM[0] = 0.;
  // CoM[1] = 0.;
  // TopGeomTools::translate(CWCForces_, CoM);
}

void DynamicPolytope::computeZMPRegion(Eigen::Vector3d comPosition)
{
  // XXX dummy zone for now: convex area formed by the polygon envelope of feet + com position
  int dim = 3;
  zmpRegion_->reset();

  // manually adding left foot points
  std::vector<Eigen::Vector3d> generators;
  auto lfPoints = robot_.surface("LeftFoot").points();
  for(auto lfPoint : lfPoints)
  {
    lfPoint = lfPoint * robot_.surface("LeftFoot").X_0_s(robot_);
    boost::shared_ptr<Generator_Rn> gn(new Generator_Rn(dim));
    boost::numeric::ublas::vector<double> coords(3);
    coords.insert_element(0, lfPoint.translation().x());
    coords.insert_element(1, lfPoint.translation().y());
    coords.insert_element(2, lfPoint.translation().z());
    gn->setCoordinates(coords);
    zmpRegion_->addGenerator(gn);
  }

  // same for right foot points
  auto rfPoints = robot_.surface("RightFoot").points();
  for(auto rfPoint : rfPoints)
  {
    rfPoint = rfPoint * robot_.surface("RightFoot").X_0_s(robot_);
    boost::shared_ptr<Generator_Rn> gn(new Generator_Rn(dim));
    boost::numeric::ublas::vector<double> coords(3);
    coords.insert_element(0, rfPoint.translation().x());
    coords.insert_element(1, rfPoint.translation().y());
    coords.insert_element(2, rfPoint.translation().z());
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

void DynamicPolytope::computeZeroMomentIntersection()
{
  zeroMomentRegion_->reset();
  // making a deep copy of the force polytope to use as base for intersection with zmp region
  // (avoids shared_ptr problems)
  politopixAPI::copyPolytope(CWCForces_, zeroMomentRegion_);

  // politopixAPI::computeIntersection(CWCForces_, zmpRegion_, zeroMomentRegion_);
  politopixAPI::computeIntersectionWithoutCheck(zeroMomentRegion_, zmpRegion_);
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

Eigen::Vector3d DynamicPolytope::computeECMP(const mc_rbdyn::Robot & robot)
{
  std::vector<std::string> contactsFSensors;
  for(auto fsensor : robot.forceSensors())
  {
    contactsFSensors.emplace_back(fsensor.name());
  }
  robotNetWrench_ = robot.netWrench(contactsFSensors);
  eCMP_ = robot.com() - (robot.com().z()) / (robot.mass() * 9.81) * robotNetWrench_.force();
  return eCMP_;
}

void DynamicPolytope::updateTrianglesGUIPolitopix()
{
  for(const auto contact : activeContacts_)
  {
    getContactMutex(forceConeTrianglesMutexes_, contact).lock();
    update3DPolyTrianglesPolitopix(frictionCones_.at(contact), forceConesTrianglesMap_.at(contact), guiScale_);
    getContactMutex(forceConeTrianglesMutexes_, contact).unlock();
    if(withMoments_)
    {
      getContactMutex(momentTrianglesMutexes_, contact).lock();
      update3DPolyTrianglesPolitopix(frictionConesMoments_.at(contact), momentPolytopesTrianglesMap_.at(contact),
                                     guiScale_);
      getContactMutex(momentTrianglesMutexes_, contact).unlock();
    }
  }
  for(const auto contact : contactsToRemove_)
  {
    getContactMutex(forceConeTrianglesMutexes_, contact).lock();
    forceConesTrianglesMap_.at(contact).clear();
    getContactMutex(forceConeTrianglesMutexes_, contact).unlock();
    if(withMoments_)
    {
      getContactMutex(momentTrianglesMutexes_, contact).lock();
      momentPolytopesTrianglesMap_.at(contact).clear();
      getContactMutex(momentTrianglesMutexes_, contact).unlock();
    }
  }

  if(!activeContacts_.empty())
  {
    // gui scale for CWC should be 1, it is position space and not force space (because eCMP)
    // update6DPolyTrianglesPolitopix(CWC_, CWCMomentTriangles_, CWCForceTriangles_, guiScale_);
    CWCForceTrianglesMutex_.lock();
    update3DPolyTrianglesPolitopix(CWCForces_, CWCForceTriangles_, 1);
    CWCForceTrianglesMutex_.unlock();
    // mc_rtc::log::info("force triangle is of size {}", CWCForceTriangles_.size());

    // scale 1 here: already position space
    ZMPTrianglesMutex_.lock();
    update3DPolyTrianglesPolitopix(zmpRegion_, ZMPTriangles_, 1);
    ZMPTrianglesMutex_.unlock();

    zeroMomentTrianglesMutex_.lock();
    update3DPolyTrianglesPolitopix(zeroMomentRegion_, zeroMomentTriangles_, 1);
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
    Normals.row(halfSpaceIndex) << halfSpace->getCoefficient(0), halfSpace->getCoefficient(1),
        halfSpace->getCoefficient(2);
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
                   mc_rtc::gui::Polyhedron(fmt::format(contact + " forces"), polyForceConfig_,
                                           [this, contact]() { return getForceConesTriangles(contact); }),
                   mc_rtc::gui::Polyhedron(fmt::format(contact + " moments"), polyMomentConfig_,
                                           [this, contact]() { return getContactMomentTriangles(contact); }));
  }

  gui.addElement(
      this, CWCCat,
      mc_rtc::gui::Polyhedron("CWC forces", polyForceConfig_, [this]() { return getCWCForceTriangles(); }),
      mc_rtc::gui::Polyhedron("CWC moments", polyMomentConfig_, [this]() { return getCWCMomentTriangles(); }),
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
