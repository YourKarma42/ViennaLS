#include <iostream>

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

#include <lsAdvect.hpp>

#include <../dev/eulerAdvect.cpp>

//#include <lsAdvect.hpp>

//____________testing end___________________________

constexpr int D = 2;
typedef double NumericType;
typedef typename lsDomain<NumericType, D>::DomainType hrleDomainType;


lsSmartPointer<lsDomain<double, D>> makeSphere(double gridDelta, double radius,
                            lsNormalizations normalization){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0.0, 0.0, 0.0};
    
    auto levelSet =
        lsSmartPointer<lsDomain<double, D>>::New(gridDelta, normalization);

    auto lsWithGeometry = lsMakeGeometry<double, D>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(origin, radius));

    lsWithGeometry.apply();


    return levelSet;

}

class velocityField : public lsVelocityField<double> {
public:
  double
  getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<double, 3>
                        & /*normalVector = hrleVectorType<double, 3>(0.)*/) {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    double velocity = 1.;
    return velocity;
  }
};

int main() {


    omp_set_num_threads(1);

    NumericType gridDelta = 0.5;

    //______________________________First____________________________________________________________________

    //std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints;

    auto start = std::chrono::high_resolution_clock::now(); 
    auto stop = std::chrono::high_resolution_clock::now(); 

    std::vector<lsSmartPointer<lsDomain<double, D>>> levelSets1;

    NumericType radius = 9.;



    //GEOMETRY CREATION
  //lsNormalizations::MANHATTEN

    lsSmartPointer<lsDomain<double, D>> levelSet1 = makeSphere(gridDelta, radius, lsNormalizations::EULER);

    levelSets1.push_back(levelSet1);  

    //ADVECTION
    auto velocities = lsSmartPointer<velocityField>::New();

    std::cout << "Advecting levelset..." << std::endl;
    eulerAdvect<NumericType, D> advectionKernel(levelSets1, velocities);

    advectionKernel.setAdvectionTime(1.);
    advectionKernel.apply();

    auto narrowband1 = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband after advection..." << std::endl;
    lsToMesh<NumericType, D>(levelSets1.back(), narrowband1, true, false, gridDelta).apply();

 
    //lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "narrowband" ).apply();
    lsVTKWriter(narrowband1, lsFileFormatEnum::VTU , "narrowband" ).apply();

    //FMM

    start = std::chrono::high_resolution_clock::now(); 

    std::cout << "Fast Marching..." << std::endl;

    //Hard coded

    //hrleVectorType<NumericType, D> origin(0.0, 0.0, 0.0);

    //lsExpandSphere<NumericType, D>(levelSets1.back(), radius, origin).apply();

    //Euler

    lsEikonalExpandTest<NumericType, D> expanderEikonal(levelSets1.back(), 5);

    expanderEikonal.apply(); 

    //Manhatten

    //lsExpand<NumericType, D>expender(levelSets1.back(), 5);

    //expender.apply();

    auto narrowband = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting FMM result..." << std::endl;
    lsToMesh<NumericType, D>(levelSets1.back(), narrowband, true, true, 4).apply();

  
    //lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "narrowband" ).apply();
    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "narrowbandFMM" ).apply();

    std::vector<std::array<NumericType, 3>> grad1;

    std::vector<NumericType> meanCurvatureSD1;

    std::vector<NumericType> meanCurvatureG;

    std::vector<NumericType> calculatedRadius;


    hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType> neighborStarIterator(levelSets1.back()->getDomain());

    hrleSparseBoxIterator<hrleDomain<NumericType, D>> neighborIterator(levelSets1.back()->getDomain(), 1);
    
    curvaturShapeDerivatives1<NumericType, D> shape1(gridDelta);

    curvaturGeneralFormula<NumericType, D> generalFormula(gridDelta);

    for(hrleConstSparseIterator<hrleDomainType> centerIt(levelSets1.back()->getDomain());
      !centerIt.isFinished(); ++centerIt){

        
        if (!centerIt.isDefined() || std::abs(centerIt.getValue()) > 0.5) {
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
          n[i] = (pos - neg) * 0.5;
          denominator += n[i] * n[i];
        }

        denominator = std::sqrt(denominator);
        if (std::abs(denominator) < 1e-12) {
          for (unsigned i = 0; i < D; ++i)
            n[i] = 0.;
        } else {
          for (unsigned i = 0; i < D; ++i) {
            //n[i] /= denominator;
          }
        }

        auto indices = centerIt.getStartIndices();

        double r = 0;

        for (unsigned i = 0; i < D; ++i) {
          r += indices[i] * gridDelta * indices[i] * gridDelta ;
        }

        calculatedRadius.push_back(std::sqrt(r) - centerIt.getValue()*gridDelta);

        meanCurvatureSD1.push_back(shape1(neighborStarIterator));

        meanCurvatureG.push_back(generalFormula(neighborIterator));

        grad1.push_back(n);       

    }

    auto narrowband2 = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(levelSets1.back(), narrowband2, true, true, 0.5).apply();
    //lsPoints



    narrowband2->insertNextVectorData(grad1, "Normals");

    narrowband2->insertNextScalarData(calculatedRadius, "radius");

    narrowband2->insertNextScalarData(meanCurvatureSD1, "curve");

    narrowband2->insertNextScalarData(meanCurvatureG, "general");

  
    lsVTKWriter(narrowband2, lsFileFormatEnum::VTU, "generalTestOutput" ).apply();





    std::cout << "Finished" << std::endl;

    return 0;
}







