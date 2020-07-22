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


    omp_set_num_threads(1);

    NumericType gridDelta = 0.25;

    //______________________________First____________________________________________________________________

    //std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints;

    std::vector<lsDomain<NumericType, D> *> levelSets;

    lsDomain<NumericType,D> levelSet = makeSphere(gridDelta, 100.);

    levelSets.push_back(&levelSet);  

    std::cout << "Converting..." << std::endl;

    lsConvertEuclid<NumericType, D>  converter(*(levelSets.back()));

    converter.apply();

    //get the active grid points of the level set
    auto activePoints = converter.getActivePoints();

    lsMesh narrowband;
    std::cout << "Extracting after conversion..." << std::endl;
    lsToMesh<NumericType, D>(levelSet, narrowband, true, true).apply(activePoints);
    //lsPoints
  
    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "ConversionOutput" ).apply();

    auto start = std::chrono::high_resolution_clock::now(); 

    std::cout << "Fast Marching..." << std::endl;

    lsEikonalExpand<NumericType, D> expander(*(levelSets.back()), activePoints);

    expander.apply(); 

    lsMesh narrowband1;
    std::cout << "Extracting after Marching..." << std::endl;
    lsToMesh<NumericType, D>(levelSet, narrowband1, true, false).apply(activePoints);
    //lsPoints
  
    lsVTKWriter(narrowband1, lsFileFormatEnum::VTU , "FMMOutput" ).apply();

    auto stop = std::chrono::high_resolution_clock::now(); 

    std::cout << "time for FMM: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << std::endl; 


    std::cout << "Finished" << std::endl;

    return 0;
}




