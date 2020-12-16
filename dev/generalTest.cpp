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

#include <../dev/lsEikonalExpandTest.hpp>



//____________testing not necessary_________________

#include <chrono>

#include <lsCalculateNormalVectors.hpp>

#include <lsPrune.hpp>

//____________testing end___________________________

constexpr int D = 3;
typedef double NumericType;
typedef typename lsDomain<NumericType, D>::DomainType hrleDomainType;


lsSmartPointer<lsDomain<double, D>> makeSphere(double gridDelta, double radius,
                            std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & lsPoints,
                            std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & narrowPoints){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0.0, 0.0, 0.0};
    
    auto levelSet =
        lsSmartPointer<lsDomain<double, D>>::New(gridDelta);

    auto lsWithGeometry = lsMakeGeometry<double, D>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(origin, radius));

    lsWithGeometry.apply();

    lsPoints = lsWithGeometry.getActivePoints();
    narrowPoints = lsWithGeometry.getNarrowPoints();


    return levelSet;

}



int main() {


    omp_set_num_threads(1);

    NumericType gridDelta = 0.5;

    //______________________________First____________________________________________________________________

    //std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints;

    auto start = std::chrono::high_resolution_clock::now(); 
    auto stop = std::chrono::high_resolution_clock::now(); 

    std::vector<lsSmartPointer<lsDomain<double, D>>> levelSets1;

    NumericType radius = 10.;


    

    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> activePoints;

    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> narrowPoints;

    

    lsSmartPointer<lsDomain<double, D>> levelSet1 = makeSphere(gridDelta, radius, activePoints, narrowPoints);

    levelSets1.push_back(levelSet1);  

    auto narrowband1 = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(levelSets1.back(), narrowband1, true, false,  gridDelta).apply();
  
    lsVTKWriter(narrowband1, lsFileFormatEnum::VTU , "/media/sf_shared/narrowbandAfterCreation" ).apply();

    start = std::chrono::high_resolution_clock::now(); 

    std::cout << "Fast Marching..." << std::endl;

    lsEikonalExpandTest<NumericType, D> expanderEikonal(levelSets1.back(), narrowPoints);

    //lsEikonalExpand<NumericType, D> expanderEikonal(levelSets1.back(), narrowPoints);

    expanderEikonal.apply(); 

    auto narrowband = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(levelSets1.back(), narrowband, true, false, gridDelta).apply();

    //lsToDiskMesh<NumericType, D>(levelSets1.back(), narrowband).apply(activePoints);
  
    //lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "narrowband" ).apply();
    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "/media/sf_shared/narrowband" ).apply();

    std::vector<std::array<NumericType, 3>> grad1;

    std::vector<NumericType> meanCurvatureSD1;

    std::vector<NumericType> meanCurvatureG;


    hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType> neighborStarIterator(levelSets1.back()->getDomain());

    hrleSparseBoxIterator<hrleDomain<NumericType, D>> neighborIterator(levelSets1.back()->getDomain(), 1);
    
    curvaturShapeDerivatives1<NumericType, D> shape1(gridDelta);

    curvaturGeneralFormula<NumericType, D> generalFormula(gridDelta);

    for(hrleConstSparseIterator<hrleDomainType> centerIt(levelSets1.back()->getDomain());
      !centerIt.isFinished(); ++centerIt){

        if (!centerIt.isDefined() || std::abs(centerIt.getValue()) > gridDelta * 0.5 ) {
        //if (!centerIt.isDefined() || std::abs(centerIt.getValue()) > (gridDelta*0.5)) {
          continue;
        } 

        // move neighborIterator to current position

        neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

        neighborStarIterator.goToIndicesSequential(centerIt.getStartIndices());

        hrleVectorType<hrleIndexType, D> tmpVec(17, 10);
        
        std::array<double, 3> n;

        NumericType denominator = 0;
        for (int i = 0; i < D; i++) {
          NumericType pos = neighborStarIterator.getNeighbor(i).getValue();
          NumericType neg = neighborStarIterator.getNeighbor(i + D).getValue();
          n[i] = (pos - neg) / (2.*gridDelta);//* 0.5;
          denominator += n[i] * n[i];
        }

        denominator = std::sqrt(denominator);
        if (std::abs(denominator) < 1e-12) {
          for (unsigned i = 0; i < D; ++i)
            n[i] = 0.;
        } else {
          for (unsigned i = 0; i < D; ++i) {
            n[i] /= denominator;
          }
        }

        meanCurvatureSD1.push_back(shape1(neighborStarIterator));

        meanCurvatureG.push_back(generalFormula(neighborIterator));

        grad1.push_back(n);       

    }

    auto narrowband2 = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(levelSets1.back(), narrowband2, true, true, gridDelta*0.5).apply();
    //lsPoints



    narrowband2->insertNextVectorData(grad1, "Normals");

    narrowband2->insertNextScalarData(meanCurvatureSD1, "curve");

    narrowband2->insertNextScalarData(meanCurvatureG, "general");

  
    lsVTKWriter(narrowband2, lsFileFormatEnum::VTU, "/media/sf_shared/normals" ).apply();





    std::cout << "Finished" << std::endl;

    return 0;
}







