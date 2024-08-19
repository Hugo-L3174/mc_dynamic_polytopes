#include "GUIComputations.h"

void update3DPolyTrianglesPolitopix(boost::shared_ptr<Polytope_Rn> & polytope,
                                    std::vector<std::array<Eigen::Vector3d, 3>> & resultTriangles,
                                    double guiScale)
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
  constIteratorOfListOfGeometricObjects<boost::shared_ptr<Generator_Rn>> generatorIter(polytope->getListOfGenerators());

  // get a point that we know is inside the polytope to compute the faces normals
  boost::numeric::ublas::vector<double> insidePoint(Rn::getDimension());

  // gravityCenter throws if there are no generators
  if(polytope->numberOfGenerators() != 0)
  {
    TopGeomTools::gravityCenter(polytope, insidePoint);
  }
  // recast it in eigen for practical reasons
  Eigen::Vector3d inside(insidePoint[0], insidePoint[1], insidePoint[2]);

  // This fills the list with vectors of the ids of the generators that compose each face (used for debug message)
  // std::vector<std::vector<unsigned int>> listOfGeneratorsPerFacet;
  // polytope->getGeneratorsPerFacet(listOfGeneratorsPerFacet);

  for(halfSpaceIter.begin(); halfSpaceIter.end() != true; halfSpaceIter.next())
  {
    // mc_rtc::log::info("Face index is {}, it has {} generators", halfSpaceIter.currentIteratorNumber(),
    //                   listOfGeneratorsPerFacet.at(halfSpaceIter.currentIteratorNumber()).size());
    std::vector<Eigen::Vector3d> vertices;
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
          Eigen::Vector3d vertex(generatorIter.current()->getCoordinate(0), generatorIter.current()->getCoordinate(1),
                                 generatorIter.current()->getCoordinate(2));
          vertices.push_back(vertex * guiScale);
          // mc_rtc::log::info("Vertex coords: {}", vertex.transpose());
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
    auto nbVertices = vertices.size();
    // mc_rtc::log::info("There are {} vertices", nbVertices);
    // XXX add nb<3 condition bc nbVerctices -2 can be negative (degen faces) and with condition on negative int is
    // evaluated to true
    for(auto i = 0; (i < nbVertices - 2) && nbVertices >= 3; i++)
    {
      // mc_rtc::log::info("Making a triangle with vertices {}, {} and {}", 0, i+1, i+2);
      // make a triangle with vertices 0, i+1, i+2 and orient and emplace it normally
      Eigen::Vector3d faceNormal;
      faceNormal = (vertices.at(i + 1) - vertices.at(0)).cross(vertices.at(i + 2) - vertices.at(0));
      faceNormal.normalize();

      auto faceOffset = halfSpaceIter.current()->getConstant();

      // testing for normal direction: if inside point of the face * normal - face offset > 0 then we need to invert
      // the face
      if(inside.dot(faceNormal) - faceOffset < 0.0)
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
