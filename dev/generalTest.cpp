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
#include <../dev/expandSphere.hpp>

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
typedef typename lsDomain<NumericType, D>::DomainType hrleDomainType;

lsDomain<double, D> makeSphere(double gridDelta, double radius){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0., 0., 0.};
    
    lsDomain<double,D> levelSet(gridDelta);

    lsMakeGeometry<double, D>(levelSet, lsSphere<double, D>(origin, radius)).apply();


    return levelSet;

}



int main() {


    omp_set_num_threads(1);

    NumericType gridDelta = 0.5;

    //______________________________First____________________________________________________________________

    //std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints;

    auto start = std::chrono::high_resolution_clock::now(); 
    auto stop = std::chrono::high_resolution_clock::now(); 

    std::vector<lsDomain<NumericType, D> *> levelSets;

    NumericType radius = 10.;

    lsDomain<NumericType,D> levelSet = makeSphere(gridDelta, radius);

    levelSets.push_back(&levelSet);  

    lsMesh pointcloud;
    std::cout << "Extracting point cloud..." << std::endl;
    lsToMesh<double, D>(*(levelSets.back()), pointcloud, true, true).apply();
    //lsToDiskMesh<double, D>(*(levelSets.back()), pointcloud).apply();
    lsVTKWriter(pointcloud, lsFileFormatEnum::VTU ,"point_cloud").apply();

/*
    int order = 1;
    lsExpand<NumericType, D>(*(levelSets.back()), 2 * (order + 2) + 1).apply();
*/


    std::cout << "Converting..." << std::endl;

    lsConvertEuclid<NumericType, D>  converter(*(levelSets.back()));

    converter.apply();

    //get the active grid points of the level set
    auto activePoints = converter.getActivePoints();

    std::cout << "active points: " << activePoints.size() << std::endl;

   // lsMesh narrowband;
    //std::cout << "Extracting after conversion..." << std::endl;
    //lsToMesh<NumericType, D>(*(levelSets.back()), narrowband, true, true).apply(activePoints);
    //lsPoints
  
    //lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "ConversionOutput" ).apply();

    start = std::chrono::high_resolution_clock::now(); 

    std::cout << "Fast Marching..." << std::endl;

    //lsEikonalExpand<NumericType, D> expander(*(levelSets.back()), activePoints);

    lsExpandSphere<NumericType, D> expander(*(levelSets.back()), activePoints, radius);


    expander.apply(); 

    std::vector<NumericType> x;
    std::vector<NumericType> y;
    std::vector<NumericType> z;


    for(hrleConstSparseIterator<hrleDomainType> it(levelSets.back()->getDomain());
        !it.isFinished(); ++it){

        x.push_back(it.getStartIndices()[0]);
        y.push_back(it.getStartIndices()[1]);

    }



    lsMesh narrowband1;
    std::cout << "Extracting after Marching..." << std::endl;
    lsToMesh<NumericType, D>(*(levelSets.back()), narrowband1, true, false).apply(activePoints);

    narrowband1.insertNextScalarData(x, "X");
    narrowband1.insertNextScalarData(y, "Y");
    //lsToMesh<NumericType, D>(levelSet, narrowband1, true, false).apply();
    //lsPoints
  
    lsVTKWriter(narrowband1, lsFileFormatEnum::VTU , "FMMOutput" ).apply();

    stop = std::chrono::high_resolution_clock::now(); 

    std::cout << "time for FMM: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << std::endl; 

    curvaturShapeDerivatives1<NumericType, D> curveCalc(gridDelta);
    std::vector<NumericType> curve;


    hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType> neighborStarIt(levelSets.back()->getDomain());

    for(hrleConstSparseIterator<hrleDomainType> it(levelSets.back()->getDomain());
        !it.isFinished(); ++it){

        //if (!it.isDefined() || std::abs(it.getValue()) > 0.5) {
        if (!it.isDefined() || (activePoints.find(it.getStartIndices()) == activePoints.end())) {
            continue;
        } 

        neighborStarIt.goToIndices(it.getStartIndices());

        curve.push_back(curveCalc(neighborStarIt));

    }



    lsMesh narrowband3;
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(*(levelSets.back()), narrowband3, true, true).apply(activePoints);
    //lsToMesh<NumericType, D>(levelSet, narrowband3, true, true).apply();
    //lsPoints

    //narrowband3.insertNextVectorData(normal, "Normal");
    narrowband3.insertNextScalarData(curve, "curvature");
  
    lsVTKWriter(narrowband3, lsFileFormatEnum::VTU , "narrowband" ).apply();


    std::cout << "Finished" << std::endl;

    return 0;
}

/*
    ITERATOR TIMING TESTS
    
    double runs = 50.;

    std::cout << "Number of runs: " << runs << std::endl;

    start = std::chrono::high_resolution_clock::now(); 

    for(int i=0; i<runs;i++){

        for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
            neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
            neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

            auto &it = neighborIt.getCenter();

            //if (!it.isDefined() || std::abs(it.getValue()) > 0.5) {
            if (!it.isDefined() || (activePoints.find(it.getStartIndices()) == activePoints.end())) {
                continue;
            } 
                  
        }
    }
    stop = std::chrono::high_resolution_clock::now(); 

    std::cout << "Star Iterator direct: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count()/runs << std::endl; 

        start = std::chrono::high_resolution_clock::now(); 

    for(int i=0; i<runs;i++){

        hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType> neighborStarIt(levelSet.getDomain());

        for(hrleConstSparseIterator<hrleDomainType> it(levelSet.getDomain());
            !it.isFinished(); ++it){

            //if (!it.isDefined() || std::abs(it.getValue()) > 0.5) {
            if (!it.isDefined() || (activePoints.find(it.getStartIndices()) == activePoints.end())) {
                continue;
            } 

            neighborStarIt.goToIndicesSequential(it.getStartIndices());

        }
    }

    stop = std::chrono::high_resolution_clock::now(); 

    std::cout << "Sparese iterator + Star Iterator goToIndicesSequential: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count()/runs << std::endl; 

        start = std::chrono::high_resolution_clock::now(); 

    for(int i=0; i<runs;i++){

        hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType> neighborStarIt(levelSet.getDomain());


        for(hrleConstSparseIterator<hrleDomainType> it(levelSet.getDomain());
            !it.isFinished(); ++it){

            //if (!it.isDefined() || std::abs(it.getValue()) > 0.5) {
            if (!it.isDefined() || (activePoints.find(it.getStartIndices()) == activePoints.end())) {
                continue;
            } 

            neighborStarIt.goToIndices(it.getStartIndices());

        }
    }

    stop = std::chrono::high_resolution_clock::now(); 

    std::cout << "Sparese iterator + Star Iterator goToIndices: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count()/runs << std::endl; 

    start = std::chrono::high_resolution_clock::now(); 

    for(int i=0; i<runs;i++){

        hrleSparseBoxIterator<hrleDomain<NumericType, D>> boxIterator(levelSet.getDomain(), 1);

        for(hrleConstSparseIterator<hrleDomainType> it(levelSet.getDomain());
            !it.isFinished(); ++it){

            //if (!it.isDefined() || std::abs(it.getValue()) > 0.5) {
            if (!it.isDefined() || (activePoints.find(it.getStartIndices()) == activePoints.end())) {
                continue;
            } 

            boxIterator.goToIndicesSequential(it.getStartIndices());

        }
    }
    stop = std::chrono::high_resolution_clock::now(); 

    std::cout << "Sparese iterator + Box Iterator goToIndicesSequential: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count()/runs << std::endl; 

        start = std::chrono::high_resolution_clock::now(); 

    for(int i=0; i<runs;i++){

        hrleSparseBoxIterator<hrleDomain<NumericType, D>> boxIterator(levelSet.getDomain(), 1);


        for(hrleConstSparseIterator<hrleDomainType> it(levelSet.getDomain());
            !it.isFinished(); ++it){

            //if (!it.isDefined() || std::abs(it.getValue()) > 0.5) {
            if (!it.isDefined() || (activePoints.find(it.getStartIndices()) == activePoints.end())) {
                continue;
            } 

            boxIterator.goToIndices(it.getStartIndices());

        }
    }
    stop = std::chrono::high_resolution_clock::now(); 

    std::cout << "Sparese iterator + Box Iterator goToIndices: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count()/runs << std::endl; 

 


*/




