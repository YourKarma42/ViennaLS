#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToMesh.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>



#include <lsEikonalExpand.hpp>
#include <../dev/expandSphere.hpp>


#include <omp.h>



#include <../dev/lsEikonalExpandTest.hpp>



//____________testing not necessary_________________

#include <chrono>

#include <lsAdvect.hpp>

#include <../dev/eulerAdvect.cpp>

//#include <lsAdvect.hpp>

//____________testing end___________________________

constexpr int D = 2;
typedef double NumericType;
typedef NumericType T;
typedef typename lsDomain<NumericType, D>::DomainType hrleDomainType;


lsSmartPointer<lsDomain<double, D>> makePlane(lsSmartPointer<lsDomain<double, D>> levelSet, double normal[3], double point[3]){



   // TODO, this is a stupid algorithm and scales with volume, which is madness
    auto &grid = levelSet->getGrid();
    hrleCoordType gridDelta = grid.getGridDelta();

    // calculate indices from sphere size
    hrleVectorType<hrleIndexType, D> index;
    hrleVectorType<hrleIndexType, D> endIndex;

    //get grid boundaries

    if(levelSet->getLevelSetNormalization() == lsNormalizations::EUCLID){

    //iterate over grid use formula

    }

    //levelSet->insertPoints(pointData);
    //levelSet->getDomain().segment();
    //levelSet->finalize(width);



    return levelSet;

}


int main() {


    omp_set_num_threads(1);

    NumericType gridDelta = 0.25;

    std::cout << "creating Plane..." << std::endl;


    std::vector<lsSmartPointer<lsDomain<double, D>>> levelSets;

    double extent = 10.;

    double bounds[2 * D] = {-extent, extent, -extent, extent};
    if (D == 3) {
        bounds[4] = -extent;
        bounds[5] = extent;
    }
    typename lsDomain<double, D>::BoundaryType boundaryCons[D];
    for (unsigned i = 0; i < D - 1; ++i) {
        boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    }
    boundaryCons[D - 1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  
    auto levelSet =
        lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta, lsNormalizations::EUCLID);

    double normal[3] = {0., 1., 0.};

    double point[3] = {0., 0., 0.};


    lsSmartPointer<lsDomain<double, D>> levelSet1 = makePlane(levelSet1, normal, point);

    levelSets.push_back(levelSet1);  


    auto narrowband1 = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband after advection..." << std::endl;
    lsToMesh<NumericType, D>(levelSets.back(), narrowband1, true, false, 1).apply();
 
    //lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "narrowband" ).apply();
    lsVTKWriter(narrowband1, lsFileFormatEnum::VTU , "narrowband" ).apply();

    //FMM


    std::cout << "Fast Marching..." << std::endl;

    lsEikonalExpandTest<NumericType, D> expanderEikonal(levelSets.back(), 5);

    expanderEikonal.apply(); 


    std::cout << "Finished" << std::endl;

    return 0;
}







