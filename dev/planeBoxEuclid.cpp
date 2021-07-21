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


typedef typename lsDomain<T, D>::PointValueVectorType pointDataType;
pointDataType pointData;


lsSmartPointer<lsDomain<double, D>> makePlane(lsSmartPointer<lsDomain<double, D>> levelSet, hrleVectorType<T, D> passedNormal, hrleVectorType<T, D> origin){



   // TODO, this is a stupid algorithm and scales with volume, which is madness
    auto &grid = levelSet->getGrid();
    hrleCoordType gridDelta = grid.getGridDelta();


    // normalise passedNormal
    double modulus = 0.;
    hrleVectorType<double, D> normal(passedNormal);
    for (unsigned i = 0; i < D; ++i) {
      modulus += normal[i] * normal[i];
    }
    modulus = std::sqrt(modulus);
    for (unsigned i = 0; i < D; ++i) {
      normal[i] /= modulus;
    }

    // check that boundary conditions are correct
    unsigned i = 0;
    bool infDimSet = false;
    for (unsigned n = 0; n < D; ++n) {
      if (grid.getBoundaryConditions(n) ==
          hrleGrid<D>::boundaryType::INFINITE_BOUNDARY) {
        if (!infDimSet) {
          i = n;
          infDimSet = true;
        } else {
          lsMessage::getInstance().addError(
              "Planes can only be created with one Infinite Boundary "
              "Condition. More than one found!");
        }
      }
    }
    if (!infDimSet) {
      lsMessage::getInstance().addError("Planes require exactly one Infinite "
                                        "Boundary Condition. None found!");
    }

    if (passedNormal[i] == 0.) {
      lsMessage::getInstance().addError(
          "lsMakeGeometry: Plane cannot be parallel to Infinite Boundary "
          "direction!");
    }


    //get grid boundaries
    //rethink width
    const T width = 2;

    if(levelSet->getLevelSetNormalization() == lsNormalizations::EUCLID){


        pointDataType pointData;

        hrleVectorType<hrleIndexType, D> index = grid.getMinIndex();
        hrleVectorType<hrleIndexType, D> minIndex = grid.getMinIndex();
        hrleVectorType<hrleIndexType, D> endIndex = grid.getMaxIndex();

        while(index<endIndex){
           
           T distToSurf = 0.;
           hrleVectorType<double, D> distPointSurf;
            
            //calculate distance

            for(int i = 0; i < D; i++){
                distPointSurf[i] = index[i] - origin[i];
            }

            for(int i = 0; i < D; i++){
                distToSurf = normal[i] * distPointSurf[i];
            }

            //add distance into the ls
            if(std::abs(distToSurf <= gridDelta)){
                pointData.push_back(std::make_pair(index, (distToSurf/gridDelta )));
                
            }

            int dim = 0;
            for (; dim < D - 1; ++dim) {
                if (index[dim] < endIndex[dim])
                break;
                index[dim] = minIndex[dim];
            }
            ++index[dim];


        }

    }

    levelSet->insertPoints(pointData);
    levelSet->getDomain().segment();
    levelSet->finalize(width);



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

    hrleVectorType<T, D> normalHrle(normal);

    hrleVectorType<T, D> pointHrle(point);


    lsSmartPointer<lsDomain<double, D>> levelSet1 = makePlane(levelSet, normalHrle, pointHrle);

    levelSets.push_back(levelSet1);  


    auto narrowband1 = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband after advection..." << std::endl;
    lsToMesh<NumericType, D>(levelSets.back(), narrowband1, true, false, 1).apply();
 
    //lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "narrowband" ).apply();
    lsVTKWriter(narrowband1, lsFileFormatEnum::VTU , "planeTEST" ).apply();

    //FMM


    std::cout << "Fast Marching..." << std::endl;

    lsEikonalExpandTest<NumericType, D> expanderEikonal(levelSets.back(), 5);

    expanderEikonal.apply(); 


    std::cout << "Finished" << std::endl;

    return 0;
}







