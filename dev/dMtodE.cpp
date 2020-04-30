#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToMesh.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>


//#include <hrleDomain.hpp>




//#include <hrleSparseBoxIterator.hpp>
//#include <hrleVectorType.hpp>

#include<lsCalculateCurvatures.hpp>



//____________testing not necessary_________________

#include <chrono>
#include <lsCalculateNormalVectors.hpp>

//____________testing end___________________________

constexpr int D = 3;
typedef double NumericType;


void createPointCloudOutput(lsDomain<NumericType, D> &passedlsDomain, std::string outputName){

  std::cout << "Calculating normal vectors..." << std::endl;




}






lsDomain<double, D> makeSphere(double gridDelta){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0., 0., 0.};
    double radius = 10.0;

    lsDomain<double,D> levelSet(gridDelta);

    lsMakeGeometry<double, D>(levelSet, lsSphere<double, D>(origin, radius)).apply();

    return levelSet;

}

constexpr double pi() { return std::atan(1)*4; }

void ManhattenToEuclidian(lsDomain<NumericType, D> &passedlsDomain){

    lsMesh mesh;

    mesh.clear();

    std::cout << "Converting Distances..." << std::endl;


    //increase lvl set size
    int order = 1;
    lsExpand<NumericType, D>(passedlsDomain, 2 * (order + 2) + 1).apply();

    std::vector<std::vector<std::array<NumericType, D>>> normalVectorsVector(
        passedlsDomain.getNumberOfSegments());

    double pointsPerSegment =
    double(2 * passedlsDomain.getDomain().getNumberOfPoints()) /
    double(passedlsDomain.getLevelSetWidth());

    auto grid = passedlsDomain.getGrid();

    //! Calculate Normalvectors
    int p = 0;
#ifdef _OPENMP
    p = omp_get_thread_num();
#endif

    std::vector<std::array<NumericType, D>> &normalVectors = normalVectorsVector[p];
    normalVectors.reserve(pointsPerSegment);

    hrleVectorType<hrleIndexType, D> startVector =
        (p == 0) ? grid.getMinGridPoint()
                : passedlsDomain.getDomain().getSegmentation()[p - 1];

    hrleVectorType<hrleIndexType, D> endVector =
      (p != static_cast<int>(passedlsDomain.getNumberOfSegments() - 1))
        ? passedlsDomain.getDomain().getSegmentation()[p]
        : grid.incrementIndices(grid.getMaxGridPoint());


    //create output vectors
    std::vector<double> oldLSValues;
    std::vector<double> newLSValues;
    

    oldLSValues.reserve(normalVectors.size());
    newLSValues.reserve(normalVectors.size());

    const NumericType gridDelta = passedlsDomain.getGrid().getGridDelta();


    //calculation of normal


    for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
            neighborIt(passedlsDomain.getDomain(), startVector);
           neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &center = neighborIt.getCenter();
        if (!center.isDefined()) {
          continue;
        } else if (std::abs(center.getValue()) > 0.5) {
          // push an empty vector to keep ordering correct
          std::array<NumericType, D> tmp = {};
          normalVectors.push_back(tmp);
          continue;
        }

        //calc normal

        //save old LS value for comparison later
        NumericType oldLSValue = center.getValue();

        std::array<NumericType, D> n;
        NumericType normN = 0.;
        for (int i = 0; i < D; i++) {
          NumericType pos = neighborIt.getNeighbor(i).getValue() - center.getValue();
          NumericType neg = center.getValue() - neighborIt.getNeighbor(i + D).getValue();
          n[i] = (pos + neg) * 0.5;
          normN += n[i] * n[i];
        }

        normN = std::sqrt(normN);

        std::array<NumericType, D> xAxis = {1. , 0. , 0.};

        //if normal points down change direction of x axis
        //if(n[0] < 0.)
            //xAxis[0] = -1.;

        //inner product
        //TODO: bad formula use different one
        NumericType ip = 0.;

        for (int i = 0; i < D; i++) {
            ip = n[i] + xAxis[i];
        }

        NumericType alphaRad = ip/normN;

        //std::cout << "|" << alphaRad*180/pi()<< "| ";

        NumericType dEuclid = center.getValue()/(std::sin(alphaRad) + std::cos(alphaRad));


        //std::cout << "|" << dEuclid << "  ";
        
        //TODO: This is tmp move somewhere else!
        //write output 

        unsigned pointId = neighborIt.getCenter().getPointId();

        // insert vertex
        std::array<unsigned, 1> vertex;
        vertex[0] = mesh.nodes.size();
        mesh.insertNextVertex(vertex);

        // insert corresponding node shifted by ls value in direction of the normal vector * LS value
        // normal vector
        std::array<double, 3> node;
        node[2] = 0.;

        //TODO: calculate Euler distance to the point we just found
        for (unsigned i = 0; i < D; ++i) {
            node[i] = double(center.getStartIndices(i)) * gridDelta;
            node[i] -= dEuclid * n[i];
        }

        mesh.insertNextNode(node);

        oldLSValues.push_back(oldLSValue);
        newLSValues.push_back(dEuclid);
    

        //normalVectors.push_back(n);
    }

    mesh.insertNextScalarData(oldLSValues, "OLD LSValues");
    mesh.insertNextScalarData(newLSValues, "NEW LSValues");

    lsVTKWriter(mesh, lsFileFormatEnum::VTU , "SpereMtoE" ).apply();


}



int main() {


    omp_set_num_threads(1);

    double gridDelta = 0.5;


    std::vector<lsDomain<double, D> *> levelSets;

    lsDomain<double,D> levelSet = makeSphere(gridDelta);

    levelSets.push_back(&levelSet);  


      

    ManhattenToEuclidian(*(levelSets.back()));


    std::cout << "Finished" << std::endl;
 

    return 0;
}



/*
lsDomain<double, D> makePlane(double gridDelta, std::vector<NumericType>& planeNormal){

      std::cout << "creating Plane..." << std::endl;

    double extent = 50;
    double bounds[2 * D] = {-extent, extent, extent, extent};
    if (D == 3) {
      bounds[4] = -extent;
      bounds[5] = extent;
    }

    typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
    for (unsigned i = 0; i < D - 1; ++i) {
      boundaryCons[i] =
          lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;

    }

    boundaryCons[D - 2] =
        lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

    lsDomain<double, D> levelSet(bounds, boundaryCons, gridDelta);


    std::vector<NumericType> origin = {0., 0., 0.};

    lsMakeGeometry<double, D>(levelSet, lsPlane<double, D>(origin, planeNormal))
      .apply();

    return levelSet;


}
*/

