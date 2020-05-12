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

lsDomain<NumericType, D> ConvertLS(lsDomain<NumericType, D> &passedlsDomain){




    std::cout << "Converting Distances..." << std::endl;


    //increase lvl set size
    int order = 1;
    lsExpand<NumericType, D>(passedlsDomain, 2 * (order + 2) + 1).apply();


    double pointsPerSegment =
    double(2 * passedlsDomain.getDomain().getNumberOfPoints()) /
    double(passedlsDomain.getLevelSetWidth());

    auto grid = passedlsDomain.getGrid();

    //! Calculate Normalvectors
    int p = 0;
#ifdef _OPENMP
    p = omp_get_thread_num();
#endif

    hrleVectorType<hrleIndexType, D> startVector =
        (p == 0) ? grid.getMinGridPoint()
                : passedlsDomain.getDomain().getSegmentation()[p - 1];

    hrleVectorType<hrleIndexType, D> endVector =
      (p != static_cast<int>(passedlsDomain.getNumberOfSegments() - 1))
        ? passedlsDomain.getDomain().getSegmentation()[p]
        : grid.incrementIndices(grid.getMaxGridPoint());

    //calculation of normal
    const NumericType gridDelta = passedlsDomain.getGrid().getGridDelta();

    //new level set with normalized values
    lsDomain<double,D> newLS(gridDelta);

    typedef typename lsDomain<NumericType, D>::PointValueVectorType pointDataType;

    pointDataType pointData;

    std::vector<double> oldVal;

    pointData.reserve(passedlsDomain.getNumberOfPoints());
    oldVal.reserve(passedlsDomain.getNumberOfPoints());

    hrleVectorType<hrleIndexType, D> index;

    for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
            neighborIt(passedlsDomain.getDomain(), startVector);
           neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &center = neighborIt.getCenter();
        if (!center.isDefined()) {
          continue;
        } else if (std::abs(center.getValue()) > 0.5) {
          // undefined run
          continue;
        }

        //save old LS value for comparison later
        NumericType oldLSValue = center.getValue();

        std::array<NumericType, D> n;
        NumericType normN = 0.;
        for (int i = 0; i < D; i++) {
          NumericType pos1 = neighborIt.getNeighbor(i).getValue();
          NumericType neg1 = neighborIt.getNeighbor(i + D).getValue();
          n[i] = (pos1 - neg1) /(2*gridDelta);
          
          normN += n[i] * n[i];
        }

        normN = std::sqrt(normN);

        //normalize vector    
        for (int i = 0; i < D; i++) {
          n[i] /= normN;
        }

        // insert corresponding node shifted by ls value in direction of the
        // normal vector

        double max = 0.;
        for (unsigned i = 0; i < D; ++i) {

            // grid position
            index[i] = center.getStartIndices(i);

            if (std::abs(n[i]) > max) {
            max = std::abs(n[i]);
            }
        }

        // the new ls value distance to the surface
        double newLsValue = center.getValue() * gridDelta * max;

        pointData.push_back(std::make_pair(index, newLsValue));

        oldVal.push_back(center.getValue());


    }
    newLS.getPointData().insertNextScalarData(oldVal, "old LS val");

    newLS.insertPoints(pointData);
    newLS.getDomain().segment();
    newLS.finalize(2);


    return newLS;
   

        
}



int main() {


    omp_set_num_threads(1);

    double gridDelta = 0.5;


    std::vector<lsDomain<NumericType, D> *> levelSets;

    lsDomain<NumericType,D> levelSet = makeSphere(gridDelta);

    levelSets.push_back(&levelSet);  

    //Converts the LS from sparse manhatten normalization to Euklid normalization
    lsDomain<NumericType,D> newLevelSet =  ConvertLS(*(levelSets.back()));

    //perform FMM

    //calculate curvatures

    //write output

    lsCalculateNormalVectors<double, D> test_normals(newLevelSet);


    lsMesh narrowband;
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(newLevelSet, narrowband).apply();

    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "test" ).apply();




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

