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
    // XXX this should be decommented but for some reason when it is, some planes are always inverted
    // The error probably comes from the sortFaceVertices function
    // hsNormal = contactPose.translation() + contactPose.rotation().transpose() * hsNormal;
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
                                    double guiScale)
{
  resultMomentTriangles.clear();
  resultForceTriangles.clear();
  // Assuming the given polytope is already computed, get the generators for each facet, then use their coordinates to
  // create the faces in mc_rtc format.

  resultMomentTriangles.reserve(polytope->numberOfHalfSpaces());
  resultForceTriangles.reserve(polytope->numberOfHalfSpaces());
  // For each half space in the polytope, get the generators that compose it
  constIteratorOfListOfGeometricObjects<boost::shared_ptr<HalfSpace_Rn>> halfSpaceIter(polytope->getListOfHalfSpaces());
  constIteratorOfListOfGeometricObjects<boost::shared_ptr<Generator_Rn>> generatorIter(polytope->getListOfGenerators());

  // get a point that we know is inside the polytope to compute the faces normals
  // Here the first 3 elements are moments, 3 after are forces: the 6d grav center makes both 3d centers
  boost::numeric::ublas::vector<double> insidePoint(Rn::getDimension());
  TopGeomTools::gravityCenter(polytope, insidePoint);
  // recast it in eigen for practical reasons
  Eigen::Vector3d insideMoment(insidePoint[0], insidePoint[1], insidePoint[2]);
  Eigen::Vector3d insideForce(insidePoint[3], insidePoint[4], insidePoint[5]);

  for(halfSpaceIter.begin(); halfSpaceIter.end() != true; halfSpaceIter.next())
  {
    // mc_rtc::log::info("Face index is {}, it has {} generators", halfSpaceIter.currentIteratorNumber(),
    //                   listOfGeneratorsPerFacet.at(halfSpaceIter.currentIteratorNumber()).size());
    std::vector<Eigen::Vector3d> verticesMoments;
    std::vector<Eigen::Vector3d> verticesForces;
    for(generatorIter.begin(); generatorIter.end() != true; generatorIter.next())
    {
      // each generator has several facets, so iterate his facets until finding the current one
      for(size_t i = 0; i < generatorIter.current()->numberOfFacets(); i++)
      {
        if(generatorIter.current()->getFacet(i) == halfSpaceIter.current())
        {
          // mc_rtc::log::info("Generator {} belongs to face {}, adding it.", generatorIter.currentIteratorNumber(),
          // halfSpaceIter.currentIteratorNumber());
          // This generator is one of the current halfspace vertex
          Eigen::Vector3d vertexMoment(generatorIter.current()->getCoordinate(0),
                                       generatorIter.current()->getCoordinate(1),
                                       generatorIter.current()->getCoordinate(2));
          Eigen::Vector3d vertexForce(generatorIter.current()->getCoordinate(3),
                                      generatorIter.current()->getCoordinate(4),
                                      generatorIter.current()->getCoordinate(5));
          verticesMoments.push_back(vertexMoment * guiScale);
          verticesForces.push_back(vertexForce * guiScale);
        }
      }
    }
    /* we got all vertices of the face in vertices vector, now order them for triangle array, ie order vertices so that
    the normal is towards the exterior
    We don't necessarily have only 3 vertices for this face! if not, more calculations are necessary to decompose
    the face into triangles
    we assume the vertices were ordered + faces are convex: then we decompose into triangles by taking the first
    vertex and making a face with the two neighbors until the second neighbor is the last vertex
    */
    auto nbVertices = verticesForces.size();
    for(auto i = 0; i < nbVertices - 2; i++)
    {
      // mc_rtc::log::info("Making a triangle with vertices {}, {} and {}", 0, i+1, i+2);
      // make a triangle with vertices 0, i+1, i+2 and orient and emplace it normally
      Eigen::Vector3d faceNormalMoments;
      Eigen::Vector3d faceNormalForces;
      faceNormalMoments =
          (verticesMoments.at(i + 1) - verticesMoments.at(0)).cross(verticesMoments.at(i + 2) - verticesMoments.at(0));
      faceNormalMoments.normalize();
      faceNormalForces =
          (verticesForces.at(i + 1) - verticesForces.at(0)).cross(verticesForces.at(i + 2) - verticesForces.at(0));
      faceNormalForces.normalize();

      auto faceOffset = halfSpaceIter.current()->getConstant();

      // testing for normal direction: if inside point of the face * normal - face offset > 0 then we need to invert
      // the face
      if(insideMoment.dot(faceNormalMoments) - faceOffset < 0.0)
      {
        resultMomentTriangles.push_back({verticesMoments.at(0), verticesMoments.at(i + 1), verticesMoments.at(i + 2)});
      }
      else
      {
        resultMomentTriangles.push_back({verticesMoments.at(0), verticesMoments.at(i + 2), verticesMoments.at(i + 1)});
      }

      if(insideForce.dot(faceNormalForces) - faceOffset < 0.0)
      {
        resultForceTriangles.push_back({verticesForces.at(0), verticesForces.at(i + 1), verticesForces.at(i + 2)});
      }
      else
      {
        resultForceTriangles.push_back({verticesForces.at(0), verticesForces.at(i + 2), verticesForces.at(i + 1)});
      }
    }
  }
}
