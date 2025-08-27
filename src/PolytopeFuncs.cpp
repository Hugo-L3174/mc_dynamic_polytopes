#include <mc_rtc/logging.h>
#include <libqhull_r/qhull_ra.h>
#include <mc_dynamic_polytopes/PolytopeFuncs.h>
#include <mc_dynamic_polytopes/WrenchCones.h>

namespace mc_dynamic_polytopes
{

void computeQhullHrep(std::vector<double> & points, boost::shared_ptr<Polytope_Rn> & polytope, int dim)
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
    HS.reset(new HalfSpace_Rn(dim));
    // Getting all coefficients of the facet normal
    // (negative for politopix convention)
    for(int coord_count = 0; coord_count < dim; coord_count++)
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

void buildFrictionConeFromContactWithHrep(int numberOfFrictionSides,
                                          const sva::PTransformd X_r1_r2,
                                          boost::shared_ptr<Polytope_Rn> & frictionCone,
                                          double m_frictionCoef,
                                          unsigned int dim)
{
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
    boost::numeric::ublas::vector<double> normal(dim);
    if(dim == 3)
    {
      normal.insert_element(0, Hrep.row(i).coeff(0));
      normal.insert_element(1, Hrep.row(i).coeff(1));
      normal.insert_element(2, Hrep.row(i).coeff(2));
    }
    if(dim == 6)
    {
      normal.insert_element(0, 0.);
      normal.insert_element(1, 0.);
      normal.insert_element(2, 0.);
      normal.insert_element(3, Hrep.row(i).coeff(0));
      normal.insert_element(4, Hrep.row(i).coeff(1));
      normal.insert_element(5, Hrep.row(i).coeff(2));
    }

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

  // XXX:(arnaud) do we need to reset the pointer
  frictionCone.reset();
  frictionCone = newCone;
}

void updatePlanesMatrixConstraint(const boost::shared_ptr<Polytope_Rn> & polytope,
                                  Eigen::MatrixXd & Normals,
                                  Eigen::VectorXd & Offsets)
{
  int nbOfPlanes = polytope->numberOfHalfSpaces();
  int dim;
  polytope->dimension() == 6 ? dim = 6 : dim = 3;
  Normals.resize(nbOfPlanes, dim);
  Offsets.resize(nbOfPlanes);

  for(int halfSpaceIndex = 0; halfSpaceIndex < nbOfPlanes; halfSpaceIndex++)
  {
    const auto halfSpace = polytope->getHalfSpace(halfSpaceIndex);
    // XXX POLITOPIX: since politopix convention seems to be inverted for Hrep inequalities
    // i.e. they check that a0 + a1*x1 + a2*x2 ... >= 0 to belong inside
    // we push inverted normals
    if(dim == 3)
    {
      Normals.row(halfSpaceIndex) << -halfSpace->getCoefficient(0), -halfSpace->getCoefficient(1),
          -halfSpace->getCoefficient(2);
    }
    else
    {
      Normals.row(halfSpaceIndex) << -halfSpace->getCoefficient(0), -halfSpace->getCoefficient(1),
          -halfSpace->getCoefficient(2), -halfSpace->getCoefficient(3), -halfSpace->getCoefficient(4),
          -halfSpace->getCoefficient(5);
    }
    Offsets(halfSpaceIndex) = halfSpace->getConstant();
  }
}

void checkAllHSInternal(const std::string & polyName, boost::shared_ptr<Polytope_Rn> & polytope)
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

} // namespace mc_dynamic_polytopes
