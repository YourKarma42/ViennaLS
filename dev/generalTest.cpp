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


#include <unordered_set>
//#include <hrleDomain.hpp>


#include <lsConvertEuclid.hpp>
#include <lsEikonalExpand.hpp>

//#include <hrleSparseBoxIterator.hpp>
//#include <hrleVectorType.hpp>



#include <lsCalculateCurvatures.hpp>

#include <omp.h>

#include <../dev/derivatives.hpp>



//____________testing not necessary_________________

#include <chrono>

#include <lsCalculateNormalVectors.hpp>

//____________testing end___________________________

constexpr int D = 2;
typedef double NumericType;

lsDomain<double, D> makeSphere(double gridDelta, double radius){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0., 0., 0.};
    
    lsDomain<double,D> levelSet(gridDelta);

    lsMakeGeometry<double, D>(levelSet, lsSphere<double, D>(origin, radius)).apply();


    return levelSet;

}



int main() {


    omp_set_num_threads(4);

    NumericType gridDelta = 0.25;

    //______________________________First____________________________________________________________________

    //std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints;

    std::vector<lsDomain<NumericType, D> *> levelSets;

    lsDomain<NumericType,D> levelSet = makeSphere(gridDelta, 10.);

    levelSets.push_back(&levelSet);  

    lsMesh pointcloud;
    std::cout << "Extracting point cloud..." << std::endl;
    //lsToMesh<double, D>(*(levelSets.back()), pointcloud).apply();
    lsToDiskMesh<double, D>(*(levelSets.back()), pointcloud).apply();
    lsVTKWriter(pointcloud, lsFileFormatEnum::VTU ,"point_cloud").apply();

/*
    std::cout << "Converting..." << std::endl;

    lsConvertEuclid<NumericType, D>  converter(*(levelSets.back()));

    converter.apply();

    //get the active grid points of the level set
    auto activePoints = converter.getActivePoints();

    std::cout << "active points: " << activePoints.size() << std::endl;

    lsMesh narrowband;
    std::cout << "Extracting after conversion..." << std::endl;
    lsToMesh<NumericType, D>(*(levelSets.back()), narrowband, true, true).apply(activePoints);
    //lsPoints
  
    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "ConversionOutput" ).apply();

*/

    auto start = std::chrono::high_resolution_clock::now(); 

    std::cout << "Fast Marching..." << std::endl;

    //lsEikonalExpand<NumericType, D> expander(*(levelSets.back()), activePoints);

    //expander.apply(); 

    int order = 1;
    lsExpand<NumericType, D>(*(levelSets.back()), 2 * (order + 2) + 1).apply();



    lsMesh narrowband1;
    std::cout << "Extracting after Marching..." << std::endl;
    //lsToMesh<NumericType, D>(levelSet, narrowband1, true, false).apply(activePoints);
    lsToMesh<NumericType, D>(levelSet, narrowband1, true, false).apply();
    //lsPoints
  
    lsVTKWriter(narrowband1, lsFileFormatEnum::VTU , "FMMOutput" ).apply();

    auto stop = std::chrono::high_resolution_clock::now(); 

    std::cout << "time for FMM: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << std::endl; 


    //curvaturShapeDerivatives1<NumericType, D> curveCalc(gridDelta);

    hrleSparseBoxIterator<hrleDomain<NumericType, D>> neighborIterator(levelSet.getDomain(), 2);


    curvaturGeneralFormula<NumericType, D> curveCalc(gridDelta);

    std::vector<NumericType> curve;



    //TODO: use different iterator
    for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
          neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
          neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

        auto &centerIt = neighborIt.getCenter();
        //if (!centerIt.isDefined() || (activePoints.find(centerIt.getStartIndices()) == activePoints.end())) {
        //  continue;
        //} 


        if (!centerIt.isDefined() || std::abs(centerIt.getValue()) > 0.5) {
          continue;
        } 

        // move neighborIterator to current position

        
        neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

        curve.push_back(curveCalc(neighborIterator));
    }


    lsMesh narrowband3;
    std::cout << "Extracting narrowband..." << std::endl;
    //lsToMesh<NumericType, D>(levelSet, narrowband3, true, true).apply(activePoints);
    lsToMesh<NumericType, D>(levelSet, narrowband3, true, true).apply();
    //lsPoints

    //narrowband.insertNextVectorData(normal, "Normal");
    narrowband3.insertNextScalarData(curve, "curvature");


  
    lsVTKWriter(narrowband3, lsFileFormatEnum::VTU , "narrowband" ).apply();


    std::cout << "Finished" << std::endl;

    return 0;
}




