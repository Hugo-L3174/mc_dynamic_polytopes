#include "GUIComputations.h"

void sortFaceVertices(std::vector<Eigen::Vector3d> & vertices, Eigen::Vector3d faceNormal)
{
  // Choose first point as reference
  Eigen::Vector3d reference = vertices[0];
  // Find a basis for the face plane to project the points in 2d on it
  Eigen::Vector3d basis1 = (vertices[1] - vertices[0]).normalized();
  Eigen::Vector3d basis2 = faceNormal.cross(basis1).normalized();

  // lambda to project on the 2d plane
  auto projectTo2D = [&basis1, &basis2, &reference](const Eigen::Vector3d & point) -> Eigen::Vector2d
  {
    Eigen::Vector3d p = point - reference;
    return Eigen::Vector2d(p.dot(basis1), p.dot(basis2));
  };
  // lambda to compare angles
  auto comparePolarAngles = [&](Eigen::Vector3d p1, Eigen::Vector3d p2, const Eigen::Vector3d reference) -> bool
  {
    Eigen::Vector2d proj1 = projectTo2D(p1);
    Eigen::Vector2d proj2 = projectTo2D(p2);

    Eigen::Vector2d refProj = projectTo2D(reference);
    Eigen::Vector2d v1 = proj1 - refProj;
    Eigen::Vector2d v2 = proj2 - refProj;

    double crossProduct = v1.x() * v2.y() - v1.y() * v2.x();

    if(crossProduct == 0)
    {
      // if colinear, sort from distance to reference
      return (v1.squaredNorm()) < (v2.squaredNorm());
    }

    // otherwise sort from cross product sign
    return crossProduct > 0;
  };

  // sort from angle order compared to the initial point
  std::sort(vertices.begin(), vertices.end(),
            [&](const Eigen::Vector3d & p1, const Eigen::Vector3d & p2)
            { return comparePolarAngles(p1, p2, reference); });
}

void update3DPolyTrianglesPolitopix(boost::shared_ptr<Polytope_Rn> & polytope,
                                    std::vector<std::array<Eigen::Vector3d, 3>> & resultTriangles,
                                    double guiScale,
                                    sva::PTransformd contactPose)
{
  // XXX CAREFUL here we fill a triangles list: this assumes we are in 3d (because each facet has dim vertices in a
  // polytope) BUT this means if we actually manipulate a 6d space the faces will be hexagons? can we assume the
  // generators are made of 2 3d matrices?
  resultTriangles.clear();
  // Assuming the given polytope is already computed, get the generators for each facet, then use their coordinates to
  // create the faces in mc_rtc format.

  resultTriangles.reserve(polytope->numberOfHalfSpaces());

  // For each half space in the polytope, get the generators that compose it
  constIteratorOfListOfGeometricObjects<boost::shared_ptr<HalfSpace_Rn>> halfSpaceIter(polytope->getListOfHalfSpaces());

  // This fills the list with vectors of the ids of the generators that compose each face
  std::vector<std::vector<unsigned int>> listOfGeneratorsPerFacet;
  polytope->getGeneratorsPerFacet(listOfGeneratorsPerFacet);

  for(halfSpaceIter.begin(); halfSpaceIter.end() != true; halfSpaceIter.next())
  {
    // get the vector of the indexes of this face's generators
    ListOfFaces listOfFacetVertices = listOfGeneratorsPerFacet.at(halfSpaceIter.currentIteratorNumber());
    // mc_rtc::log::info("Face index is {}, it has {} generators", halfSpaceIter.currentIteratorNumber(),
    //                   listOfFacetVertices.size());

    std::vector<Eigen::Vector3d> vertices;

    // get vertices coordinates of this face to build the triangles composing it
    for(int faceGens = 0; faceGens < listOfFacetVertices.size(); faceGens++)
    {
      int generatorIndex = listOfFacetVertices.at(faceGens);
      // mc_rtc::log::info("Generator {} belongs to face {}, adding it.", generatorIndex,
      //                   halfSpaceIter.currentIteratorNumber());
      // This generator is one of the current halfspace vertex
      Eigen::Vector3d vertex(polytope->getGenerator(generatorIndex)->getCoordinate(0),
                             polytope->getGenerator(generatorIndex)->getCoordinate(1),
                             polytope->getGenerator(generatorIndex)->getCoordinate(2));
      vertex *= guiScale;
      // vertex pose is X_c_v (contact frame) but we want it in world frame
      // X_0_v = X_c_v * X_0_c
      // 3d pose rotated will be pt.r_ + pt.E_.transpose() * r_
      // vertex = contactPose.translation() + contactPose.rotation().transpose() * vertex;
      vertex = (sva::PTransformd(vertex) * contactPose).translation();
      vertices.push_back(vertex);

      // mc_rtc::log::info("Vertex coords: {}", vertex.transpose());
    }

    /* we got all vertices of the face in vertices vector, now order them for triangle array, ie order vertices so that
    the normal is towards the exterior
    We don't necessarily have only 3 vertices for this face! if not, more calculations are necessary to decompose
    the face into triangles
    we assume the vertices were ordered + faces are convex: then we decompose into triangles by taking the first
    vertex and making a face with the two neighbors until the second neighbor is the last vertex
    */
    auto nbVertices = vertices.size();
    // mc_rtc::log::info("There are {} vertices", nbVertices);

    // XXX POLITOPIX: since politopix convention seems to be inverted for Hrep inequalities
    // i.e. they check that a0 + a1*x1 + a2*x2 ... >= 0 to belong inside
    // we push inverted normals
    Eigen::Vector3d hsNormal(-halfSpaceIter.current()->getCoefficient(0), -halfSpaceIter.current()->getCoefficient(1),
                             -halfSpaceIter.current()->getCoefficient(2));
    // rotate hsNormal to world frame as well
    hsNormal = contactPose.rotation().transpose() * hsNormal;
    hsNormal.normalize();

    // Safety check: sometimes facets have zero vertices idk why, probably degenerated faces. The sorting throws if this
    // is the case
    if(nbVertices != 0)
    {
      // ordering the vertices
      sortFaceVertices(vertices, hsNormal);
    }

    // add nbVertices<3 condition because nbVertices -2 can be negative (degen faces) and with condition on negative int
    // is evaluated to true
    for(auto i = 0; (i < nbVertices - 2) && (nbVertices >= 3); i++)
    {
      // mc_rtc::log::info("Making a triangle with vertices {}, {} and {}", 0, i+1, i+2);
      // make a triangle with vertices 0, i+1, i+2 and orient and emplace it normally
      Eigen::Vector3d faceNormal;
      faceNormal = (vertices.at(i + 1) - vertices.at(0)).cross(vertices.at(i + 2) - vertices.at(0));
      // faceNormal *= -1.0;
      faceNormal.normalize();
      // mc_rtc::log::info("computed normal is {}", faceNormal);
      // mc_rtc::log::info("hsNormal is {}", hsNormal);

      // testing for normal direction: if normal of the triangle face * normal of the facet < 0 then we need to invert
      // the face (politopix and gui conventions are inverted I think)
      if(faceNormal.dot(hsNormal) > 0.0)
      {
        resultTriangles.push_back({vertices.at(0), vertices.at(i + 1), vertices.at(i + 2)});
      }
      else
      {
        resultTriangles.push_back({vertices.at(0), vertices.at(i + 2), vertices.at(i + 1)});
      }
    }
  }
}

void update6DPolyTrianglesPolitopix(boost::shared_ptr<Polytope_Rn> & polytope,
                                    std::vector<std::array<Eigen::Vector3d, 3>> & resultMomentTriangles,
                                    std::vector<std::array<Eigen::Vector3d, 3>> & resultForceTriangles,
                                    double guiScale,
                                    sva::PTransformd contactPose)
{
  resultMomentTriangles.clear();
  resultForceTriangles.clear();
  // Assuming the given polytope is already computed, get the generators for each facet, then use their coordinates to
  // create the faces in mc_rtc format.

  resultMomentTriangles.reserve(polytope->numberOfHalfSpaces());
  resultForceTriangles.reserve(polytope->numberOfHalfSpaces());
  // For each half space in the polytope, get the generators that compose it
  constIteratorOfListOfGeometricObjects<boost::shared_ptr<HalfSpace_Rn>> halfSpaceIter(polytope->getListOfHalfSpaces());

  // This fills the list with vectors of the ids of the generators that compose each face
  std::vector<std::vector<unsigned int>> listOfGeneratorsPerFacet;
  polytope->getGeneratorsPerFacet(listOfGeneratorsPerFacet);

  for(halfSpaceIter.begin(); halfSpaceIter.end() != true; halfSpaceIter.next())
  {
    // get the vector of the indexes of this face's generators
    ListOfFaces listOfFacetVertices = listOfGeneratorsPerFacet.at(halfSpaceIter.currentIteratorNumber());
    // mc_rtc::log::info("Face index is {}, it has {} generators", halfSpaceIter.currentIteratorNumber(),
    //                   listOfFacetVertices.size());

    std::vector<Eigen::Vector3d> forceVertices;
    std::vector<Eigen::Vector3d> momentVertices;

    // get vertices coordinates of this face to build the triangles composing it
    for(int faceGens = 0; faceGens < listOfFacetVertices.size(); faceGens++)
    {
      int generatorIndex = listOfFacetVertices.at(faceGens);
      // mc_rtc::log::info("Generator {} belongs to face {}, adding it.", generatorIndex,
      //                   halfSpaceIter.currentIteratorNumber());
      // This generator is one of the current halfspace vertex
      Eigen::Vector3d momentVertex(polytope->getGenerator(generatorIndex)->getCoordinate(0),
                                   polytope->getGenerator(generatorIndex)->getCoordinate(1),
                                   polytope->getGenerator(generatorIndex)->getCoordinate(2));

      Eigen::Vector3d forceVertex(polytope->getGenerator(generatorIndex)->getCoordinate(3),
                                  polytope->getGenerator(generatorIndex)->getCoordinate(4),
                                  polytope->getGenerator(generatorIndex)->getCoordinate(5));

      momentVertex *= guiScale;
      forceVertex *= guiScale;

      // vertex pose is X_c_v (contact frame) but we want it in world frame
      // X_0_v = X_c_v * X_0_c
      // 3d pose rotated will be pt.r_ + pt.E_.transpose() * r_
      // vertex = contactPose.translation() + contactPose.rotation().transpose() * vertex;
      momentVertex = (sva::PTransformd(momentVertex) * contactPose).translation();
      momentVertices.push_back(momentVertex);

      forceVertex = (sva::PTransformd(forceVertex) * contactPose).translation();
      forceVertices.push_back(momentVertex);

      // mc_rtc::log::info("Vertex coords: {}", vertex.transpose());
    }

    /* we got all vertices of the face in vertices vector, now order them for triangle array, ie order vertices so that
    the normal is towards the exterior
    We don't necessarily have only 3 vertices for this face! if not, more calculations are necessary to decompose
    the face into triangles
    we assume the vertices were ordered + faces are convex: then we decompose into triangles by taking the first
    vertex and making a face with the two neighbors until the second neighbor is the last vertex
    */
    auto nbForceVertices = forceVertices.size();
    auto nbMomentVertices = momentVertices.size();
    // mc_rtc::log::info("There are {} vertices", nbVertices);

    // XXX POLITOPIX: since politopix convention seems to be inverted for Hrep inequalities
    // i.e. they check that a0 + a1*x1 + a2*x2 ... >= 0 to belong inside
    // we push inverted normals
    Eigen::Vector3d hsNormalMoment(-halfSpaceIter.current()->getCoefficient(0),
                                   -halfSpaceIter.current()->getCoefficient(1),
                                   -halfSpaceIter.current()->getCoefficient(2));

    Eigen::Vector3d hsNormalForce(-halfSpaceIter.current()->getCoefficient(3),
                                  -halfSpaceIter.current()->getCoefficient(4),
                                  -halfSpaceIter.current()->getCoefficient(5));
    // rotate hsNormal to world frame as well
    hsNormalMoment = contactPose.rotation().transpose() * hsNormalMoment;
    hsNormalMoment.normalize();
    hsNormalForce = contactPose.rotation().transpose() * hsNormalForce;
    hsNormalForce.normalize();

    // Safety check: sometimes facets have zero vertices idk why, probably degenerated faces. The sorting throws if this
    // is the case
    if(nbMomentVertices != 0)
    {
      // ordering the vertices
      sortFaceVertices(momentVertices, hsNormalMoment);
    }
    if(nbForceVertices != 0)
    {
      // ordering the vertices
      sortFaceVertices(forceVertices, hsNormalForce);
    }

    // add nbVertices<3 condition because nbVertices -2 can be negative (degen faces) and with condition on negative int
    // is evaluated to true
    for(auto i = 0; (i < nbMomentVertices - 2) && (nbMomentVertices >= 3); i++)
    {
      // mc_rtc::log::info("Making a triangle with vertices {}, {} and {}", 0, i+1, i+2);
      // make a triangle with vertices 0, i+1, i+2 and orient and emplace it normally
      Eigen::Vector3d faceNormal;
      faceNormal =
          (momentVertices.at(i + 1) - momentVertices.at(0)).cross(momentVertices.at(i + 2) - momentVertices.at(0));
      // faceNormal *= -1.0;
      faceNormal.normalize();
      // mc_rtc::log::info("computed normal is {}", faceNormal);
      // mc_rtc::log::info("hsNormal is {}", hsNormal);

      // testing for normal direction: if normal of the triangle face * normal of the facet < 0 then we need to invert
      // the face (politopix and gui conventions are inverted I think)
      if(faceNormal.dot(hsNormalMoment) > 0.0)
      {
        resultMomentTriangles.push_back({momentVertices.at(0), momentVertices.at(i + 1), momentVertices.at(i + 2)});
      }
      else
      {
        resultMomentTriangles.push_back({momentVertices.at(0), momentVertices.at(i + 2), momentVertices.at(i + 1)});
      }
    }

    for(auto i = 0; (i < nbForceVertices - 2) && (nbForceVertices >= 3); i++)
    {
      // mc_rtc::log::info("Making a triangle with vertices {}, {} and {}", 0, i+1, i+2);
      // make a triangle with vertices 0, i+1, i+2 and orient and emplace it normally
      Eigen::Vector3d faceNormal;
      faceNormal = (forceVertices.at(i + 1) - forceVertices.at(0)).cross(forceVertices.at(i + 2) - forceVertices.at(0));
      // faceNormal *= -1.0;
      faceNormal.normalize();
      // mc_rtc::log::info("computed normal is {}", faceNormal);
      // mc_rtc::log::info("hsNormal is {}", hsNormal);

      // testing for normal direction: if normal of the triangle face * normal of the facet < 0 then we need to invert
      // the face (politopix and gui conventions are inverted I think)
      if(faceNormal.dot(hsNormalForce) > 0.0)
      {
        resultForceTriangles.push_back({forceVertices.at(0), forceVertices.at(i + 1), forceVertices.at(i + 2)});
      }
      else
      {
        resultForceTriangles.push_back({forceVertices.at(0), forceVertices.at(i + 2), forceVertices.at(i + 1)});
      }
    }
  }
}
