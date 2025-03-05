#include "WrenchCones.h"

// Julien's method, not convinced
Eigen::MatrixXd linearizedFrictionCone(int numberOfFrictionSides, Eigen::Matrix3d m_rotation, double m_frictionCoef)
{
  Eigen::MatrixXd F(numberOfFrictionSides + 2, 3); // XXX something wrong with the +2 ? why ?
  Eigen::RowVector3d line;

  //   XXX why are the first 2 lines of the cone matrix 0 0 1 * rot^T and 0 0 -1 * rot^T
  line << 0, 0, 1;
  line = line * (m_rotation.transpose());
  F.row(0) = line;

  line << 0, 0, -1;
  line = line * (m_rotation.transpose());
  F.row(1) = line;

  double Dtheta(2 * M_PI / numberOfFrictionSides);
  // std::cout << "Dtheta: " << Dtheta << '\n';
  for(int i = 0; i < numberOfFrictionSides; ++i)
  {
    line << cos(i * Dtheta), sin(i * Dtheta), -m_frictionCoef * cos(Dtheta / 2);
    // line << cos(i*Dtheta), sin(i*Dtheta), -m_frictionCoef;
    line = line * (m_rotation.transpose());
    F.row(i + 2) = line;
  }

  return F;
}

std::vector<Eigen::Vector3d> generatePolyhedralConeGens(int numberOfFrictionSides,
                                                        Eigen::Matrix3d m_rotation,
                                                        double m_frictionCoef,
                                                        double maxForce)
{
  std::vector<Eigen::Vector3d> generators(numberOfFrictionSides);
  // here unit vector, but can be capped for cone?
  Eigen::Vector3d normal(Eigen::Vector3d::UnitZ());
  Eigen::Vector3d tan(Eigen::Vector3d::UnitX());
  double angle = std::atan(m_frictionCoef);

  // gen is the max tangential axis tolerated around the normal of the contact, deduced from the friction coeff
  Eigen::Vector3d gen = Eigen::AngleAxisd(angle, tan) * normal;
  // XXX as we cheat for now and build a polytope instead of a polyhedral cone, we scale the generating vector to have a
  // normal value of maxForce
  gen *= (maxForce / gen.z());
  // step is the scale decomposition (precision) with which to compute the actual cone (linearization)
  double step = (M_PI * 2.) / numberOfFrictionSides;

  for(int i = 0; i < numberOfFrictionSides; i++)
  {
    // each generator is formed by the limit points of the linearized cone around the contact normal
    generators[i] = m_rotation.transpose() * Eigen::AngleAxisd(step * i, normal) * gen;
    // mc_rtc::log::info("generator {} of this cone is {}", i, generators[i].transpose());
  }
  return generators;
}

Eigen::MatrixXd generatePolyhedralConeHRep(int numberOfFrictionSides, Eigen::Matrix3d rotX_r1_r2, double m_frictionCoef)
{
  Eigen::MatrixXd HRep(numberOfFrictionSides, 3);
  Eigen::Vector3d contactNormal(Eigen::Vector3d::UnitZ());
  Eigen::Vector3d tan(Eigen::Vector3d::UnitX());

  // The angle to the contact normal of the friction cone is atan(mu)
  // (mu for external approximation, mu/sqrt(2) for internal, let's pick internal for H-rep)
  // But here for hrep we want the normals of the linearized cone's faces
  // --> there is a 90° angle to add to get the face normal
  double angle = (M_PI / 2.) + atan(m_frictionCoef / sqrt(2));
  // This is the first face normal
  Eigen::Vector3d normal = Eigen::AngleAxisd(angle, tan) * contactNormal;

  // step is the scale decomposition (precision) with which to compute the actual cone (linearization)
  double step = (M_PI * 2.) / numberOfFrictionSides;

  // here we compute the hrep: the rows will be the normals of the cone faces in the controlled frame
  for(int i = 0; i < numberOfFrictionSides; i++)
  {
    // Rotation around contact normal for each decomposed angle, then transposed in the relative contact frame
    HRep.row(i) = rotX_r1_r2.transpose() * Eigen::AngleAxisd(step * i, contactNormal) * normal;
  }
  return HRep;
}

Eigen::MatrixXd compute6DGeneratorsMatrixSingleCone(
    Eigen::Vector3d applicationPoint,
    int numberOfFrictionSides,
    std::pair<std::pair<double, double>, sva::PTransformd> contactSurface,
    double m_frictionCoef)
{
  Eigen::MatrixXd genMatrix;
  genMatrix.resize(6, numberOfFrictionSides
                          * 4); // 4 is fixed because we assume rectangular contacts (4 individual contact points)
  Eigen::Index col(0);

  // mc_rtc::log::info("generating cones for a contact point");
  auto generators =
      generatePolyhedralConeGens(numberOfFrictionSides, contactSurface.second.rotation(), m_frictionCoef, 1);
  // contactPoint first is the pair of xHalfLength and yHalfLength of the rectangular contact
  std::vector<Eigen::Vector3d> points;
  points.emplace_back(contactSurface.first.first, contactSurface.first.second, 0);
  points.emplace_back(contactSurface.first.first, -contactSurface.first.second, 0);
  points.emplace_back(-contactSurface.first.first, -contactSurface.first.second, 0);
  points.emplace_back(-contactSurface.first.first, contactSurface.first.second, 0);
  // mc_rtc::log::info("computing matrix using cone generators:");
  for(auto p : points)
  {
    // XXX r is the extremity point of the surface : center + offset using half dimensions (-application point in world
    // to get desired frame)
    Eigen::Vector3d r = contactSurface.second.translation() + p - applicationPoint;
    // mc_rtc::log::info("point half dim offset is {}, contact point associated is {}", p.transpose(), r.transpose());
    for(auto g : generators)
    {
      // here we stack the matrix using the "limits" of the surface contact (points r) and associate each generator
      // found for the contact this is the angular part of the 6d vector (resulting moment of the force generator at
      // application point)
      genMatrix.col(col).segment<3>(0).noalias() = skewMatrix(r) * g;
      // this is the translational part
      genMatrix.col(col).segment<3>(3) = g;
      // mc_rtc::log::info("applying generator to contact point : 6d vect is {}", genMatrix.col(col).transpose());
      col += 1;
    }
  }
  return genMatrix;
}

Eigen::MatrixXd compute6DGeneratorsMatrixRaysCones(
    Eigen::Vector3d applicationPoint,
    int numberOfFrictionSides,
    std::vector<std::pair<std::pair<double, double>, sva::PTransformd>> contactSurfaces,
    double m_frictionCoef)
{
  Eigen::MatrixXd genMatrix;
  genMatrix.resize(6, numberOfFrictionSides
                          * 4); // 4 is fixed because we assume rectangular contacts (4 individual contact points)
  Eigen::Index col(0);

  for(auto contactSurface : contactSurfaces)
  {
    // mc_rtc::log::info("generating cones for a contact point");
    auto generators =
        generatePolyhedralConeGens(numberOfFrictionSides, contactSurface.second.rotation(), m_frictionCoef, 1);
    // contactPoint first is the pair of xHalfLength and yHalfLength of the rectangular contact
    std::vector<Eigen::Vector3d> points;
    points.emplace_back(contactSurface.first.first, contactSurface.first.second, 0);
    points.emplace_back(contactSurface.first.first, -contactSurface.first.second, 0);
    points.emplace_back(-contactSurface.first.first, -contactSurface.first.second, 0);
    points.emplace_back(-contactSurface.first.first, contactSurface.first.second, 0);
    // mc_rtc::log::info("computing matrix using cone generators:");
    for(auto p : points)
    {
      // XXX r is the extremity point of the surface : center + offset using half dimensions (-application point in
      // world to get desired frame)
      Eigen::Vector3d r = contactSurface.second.translation() + p - applicationPoint;
      // mc_rtc::log::info("point half dim offset is {}, contact point associated is {}", p.transpose(), r.transpose());
      for(auto g : generators)
      {
        // here we stack the matrix using the "limits" of the surface contact (points r) and associate each generator
        // found for the contact this is the angular part of the 6d vector (resulting moment of the force generator at
        // application point)
        genMatrix.col(col).segment<3>(0).noalias() = skewMatrix(r) * g;
        // this is the translational part
        genMatrix.col(col).segment<3>(3) = g;
        // mc_rtc::log::info("applying generator to contact point : 6d vect is {}", genMatrix.col(col).transpose());
        col += 1;
      }
    }
  }
  return genMatrix;
}

Eigen::Matrix3d skewMatrix(const Eigen::Vector3d v)
{
  Eigen::Matrix3d mat;
  mat <<
      // clang-format off
        0.,  -v(2), v(1),
        v(2),   0.,-v(0),
       -v(1), v(0),   0.;
  // clang-format on
  return mat;
}

void findHalfWidthLength(const mc_rbdyn::Surface & surface, double & halfWidth, double & halfLength)
{
  const auto & surfacePoints = surface.points();
  // Find boundaries in surface frame along the surface's sagital (x) and lateral (y) direction
  double minSagital = std::numeric_limits<double>::max();
  double minLateral = std::numeric_limits<double>::max();
  double maxSagital = -std::numeric_limits<double>::max();
  double maxLateral = -std::numeric_limits<double>::max();
  for(const auto & point : surfacePoints)
  {
    // Points are defined in body frame, convert to surface frame
    Eigen::Vector3d surfacePoint = surface.X_b_s().rotation() * (point.translation() - surface.X_b_s().translation());
    double x = surfacePoint.x();
    double y = surfacePoint.y();
    minSagital = std::min(minSagital, x);
    maxSagital = std::max(maxSagital, x);
    minLateral = std::min(minLateral, y);
    maxLateral = std::max(maxLateral, y);
  }

  halfLength = (maxSagital - minSagital) / 2.;
  halfWidth = (maxLateral - minLateral) / 2.;
}
