
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

#include <../dev/expandSphere.hpp>

#include <../dev/eulerAdvect.cpp>



//____________testing not necessary_________________

#include <chrono>

#include <lsCalculateNormalVectors.hpp>

constexpr int D = 3 ;
typedef double NumericType;

typedef typename lsDomain<NumericType, D>::DomainType hrleDomainType;



lsSmartPointer<lsDomain<double, D>> makeSphere(double gridDelta, double radius,
                            std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & narrowPoints){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0.0, 0.0, 0.0};
    
    auto levelSet =
        lsSmartPointer<lsDomain<double, D>>::New(gridDelta);

    auto lsWithGeometry = lsMakeGeometry<double, D>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(origin, radius));

    lsWithGeometry.apply();

    narrowPoints = lsWithGeometry.getNarrowPoints();


    return levelSet;

}

void create_output_Manhatten(lsSmartPointer<lsDomain<double, D>> levelSet,
  std::string output_name){

    hrleSparseBoxIterator<hrleDomain<NumericType, D>> neighborIterator(levelSet->getDomain(), 1);

    std::vector<std::array<NumericType, 3>> normal;

    std::vector<std::array<NumericType, 3>> secondOrderDerivatives1;

    std::vector<std::array<NumericType, 3>> secondOrderDerivatives2;



    std::vector<NumericType> grad1;

    std::vector<NumericType> meanCurvatureGeneralFormula;

    std::vector<NumericType> meanCurvatureGeneralFormulaBig;

    std::vector<NumericType> meanCurvatureGeneralFormulaBigBias;

    std::vector<NumericType> meanCurvatureSD1;

    std::vector<NumericType> grad2;

    std::vector<NumericType> meanCurvature2;

    std::vector<NumericType> meanCurvatureSD2;

    std::vector<NumericType> meanCurvature3;

    std::vector<NumericType> meanCurvatureNew;

    std::vector<NumericType> curveTest;

    std::vector<NumericType> meanCurveShapeBias;


    //std::vector<std::array<NumericType, D>> normal;


    NumericType gridDelta = levelSet->getGrid().getGridDelta();

    variationOfNormals<NumericType, D> variationOfGrad(gridDelta);

    curvaturShapeDerivatives1<NumericType, D> shape1(gridDelta);

    curvaturShapeDerivatives2<NumericType, D> shape2(gridDelta);

    curvaturGeneralFormula<NumericType, D> generalFormula(gridDelta);

    curvaturShapeBias<NumericType, D> shapeBias(gridDelta);

    curvaturGeneralFormulaBigStencil <NumericType, D> generalFormulaBig(gridDelta);

    curvaturGeneralFormulaBigStencilBias <NumericType, D> generalFormulaBigBias(gridDelta);

    curvaturTest<NumericType, D> test(gridDelta);


  
   //     ? passedlsDomain.getDomain().getSegmentation()[p]
   //     : grid.incrementIndices(grid.getMaxGridPoint());


    //TODO: use hrleConstSparseIterator<hrleDomainType>
    for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
          neighborIt(levelSet->getDomain(), levelSet->getGrid().getMinGridPoint());
          neighborIt.getIndices() < levelSet->getGrid().incrementIndices(levelSet->getGrid().getMaxGridPoint()); neighborIt.next()) {

        auto &centerIt = neighborIt.getCenter();
        if (!centerIt.isDefined() || (std::abs(centerIt.getValue()) > 0.5)) {
          continue;
        } 

        // move neighborIterator to current position
        neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

        meanCurvatureSD1.push_back(shape1(neighborIt));

        meanCurvatureSD2.push_back(shape2(neighborIterator));

        meanCurveShapeBias.push_back(shapeBias(neighborIterator));

        meanCurvatureNew.push_back(variationOfGrad(neighborIterator));

        meanCurvatureGeneralFormula.push_back(generalFormula(neighborIterator));

        meanCurvatureGeneralFormulaBig.push_back(generalFormulaBig(neighborIterator));

        meanCurvatureGeneralFormulaBigBias.push_back(generalFormulaBigBias(neighborIterator));

        curveTest.push_back(test(neighborIterator));

        

    }

    auto narrowband = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(levelSet, narrowband, true, true).apply();
    //lsPoints


    //narrowband.insertNextVectorData(normal, "Normal");
    narrowband->insertNextScalarData(meanCurvatureSD1, "shape operator small stencil");
    narrowband->insertNextScalarData(meanCurvatureSD2, "shape operator big stencil");
    narrowband->insertNextScalarData(meanCurvatureGeneralFormula, "general formula");
    narrowband->insertNextScalarData(meanCurvatureGeneralFormulaBig, "general formula big stencil");
    narrowband->insertNextScalarData(meanCurvatureGeneralFormulaBigBias, "general formula big stencil bias");
    narrowband->insertNextScalarData(curveTest, "Curvature Test");
    narrowband->insertNextScalarData(meanCurvatureNew, "variation of normals");
    narrowband->insertNextScalarData(meanCurveShapeBias, "shape operator Bias");

  
    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , output_name ).apply();


  }

void create_output_Euklid(lsSmartPointer<lsDomain<double, D>> levelSet,
  std::string output_name){

    //std::vector<NumericType> distance;


    std::vector<NumericType> grad1;

    std::vector<NumericType> meanCurvatureGeneralFormula;

    std::vector<NumericType> meanCurvatureGeneralFormulaBig;

    std::vector<NumericType> meanCurvatureGeneralFormulaBigBias;

    std::vector<NumericType> meanCurvatureSD1;

    std::vector<NumericType> grad2;

    std::vector<NumericType> meanCurvature2;

    std::vector<NumericType> meanCurvatureSD2;

    std::vector<NumericType> meanCurvature3;

    std::vector<NumericType> meanCurvatureNew;

    std::vector<NumericType> curveTest;

    std::vector<NumericType> meanCurveShapeBias;


    //std::vector<std::array<NumericType, D>> normal;


    NumericType gridDelta = levelSet->getGrid().getGridDelta();


    lsEikonalExpandTest<NumericType, D> expanderEikonal(levelSet, 5);

    expanderEikonal.apply(); 


    variationOfNormals<NumericType, D> variationOfGrad(gridDelta);

    curvaturShapeDerivatives1<NumericType, D> shape1(gridDelta);

    curvaturShapeDerivatives2<NumericType, D> shape2(gridDelta);

    curvaturGeneralFormula<NumericType, D> generalFormula(gridDelta);

    curvaturShapeBias<NumericType, D> shapeBias(gridDelta);

    curvaturGeneralFormulaBigStencil <NumericType, D> generalFormulaBig(gridDelta);

    curvaturGeneralFormulaBigStencilBias <NumericType, D> generalFormulaBigBias(gridDelta);

    curvaturTest<NumericType, D> test(gridDelta);


  
   //     ? passedlsDomain.getDomain().getSegmentation()[p]
   //     : grid.incrementIndices(grid.getMaxGridPoint());

    hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType> neighborStarIterator(levelSet->getDomain());

    hrleSparseBoxIterator<hrleDomain<NumericType, D>> neighborIterator(levelSet->getDomain(), 1);


    for(hrleConstSparseIterator<hrleDomainType> centerIt(levelSet->getDomain());
      !centerIt.isFinished(); ++centerIt){
   // for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
    //      neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
     //     neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

        //auto &centerIt = neighborIt.getCenter();
        //if (!centerIt.isDefined() || (lsPoints.find(centerIt.getStartIndices()) == lsPoints.end())) {
        if (!centerIt.isDefined() || std::abs(centerIt.getValue()) > gridDelta* 0.5) {
          continue;
        } 

        // move neighborIterator to current position

        neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

        neighborStarIterator.goToIndicesSequential(centerIt.getStartIndices());
        

        meanCurvatureSD1.push_back(shape1(neighborStarIterator));

        meanCurvatureSD2.push_back(shape2(neighborIterator));

        meanCurveShapeBias.push_back(shapeBias(neighborIterator));

        meanCurvatureNew.push_back(variationOfGrad(neighborIterator));

        meanCurvatureGeneralFormula.push_back(generalFormula(neighborIterator));

        meanCurvatureGeneralFormulaBig.push_back(generalFormulaBig(neighborIterator));

        meanCurvatureGeneralFormulaBigBias.push_back(generalFormulaBigBias(neighborIterator));

        curveTest.push_back(test(neighborIterator));

        

    }

    auto narrowband = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(levelSet, narrowband, true, true, gridDelta*0.5).apply();
    //lsPoints


    //narrowband.insertNextVectorData(normal, "Normal");
    narrowband->insertNextScalarData(meanCurvatureSD1, "shape operator small stencil");
    narrowband->insertNextScalarData(meanCurvatureSD2, "shape operator big stencil");
    narrowband->insertNextScalarData(meanCurvatureGeneralFormula, "general formula");
    narrowband->insertNextScalarData(meanCurvatureGeneralFormulaBig, "general formula big stencil");
    narrowband->insertNextScalarData(meanCurvatureGeneralFormulaBigBias, "general formula big stencil bias");
    narrowband->insertNextScalarData(curveTest, "Curvature Test");
    narrowband->insertNextScalarData(meanCurvatureNew, "variation of normals");
    narrowband->insertNextScalarData(meanCurveShapeBias, "shape operator Bias");

  
    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , output_name ).apply();


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



int main(int argc, char* argv[]) {

    int numThreads = 1;

    NumericType gridDelta = 0.5;

    if(argc != 1){
        numThreads = std::stoi(argv[1]);
    }

    std::cout << "running program with " << numThreads << " threads" << std::endl;


    omp_set_num_threads(numThreads);

   std::vector<lsSmartPointer<lsDomain<double, D>>> levelSets1;

    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> narrowPoints;


    NumericType radius = 95.;

    lsSmartPointer<lsDomain<double, D>> levelSet = makeSphere(gridDelta, radius, narrowPoints);

    levelSets1.push_back(levelSet);  

    auto velocities = lsSmartPointer<velocityField>::New();

    std::cout << "Advecting levelset..." << std::endl;
    eulerAdvect<NumericType, D> advectionKernel(levelSets1, velocities);


    advectionKernel.setAdvectionTime(5.);
    advectionKernel.apply();

//    auto mesh = lsSmartPointer<lsMesh>::New();
//    lsToDiskMesh<double, D>(levelSet, mesh).apply();
//    lsVTKWriter(mesh, "sphere-" + std::to_string(1) + ".vtk").apply();

    std::cout << "Calculating Curvatures ..." << std::endl;

    create_output_Euklid(levelSets1.back(), "advectedEuler");
    //create_output_Manhatten(levelSet,"advectedManhatten");

/*
    unsigned counter = 1;
    for (double time = 0; time < 1.; time += advectionKernel.getAdvectedTime()) {
        advectionKernel.apply();

        auto mesh = lsSmartPointer<lsMesh>::New();
        lsToDiskMesh<double, D>(newLayer, mesh).apply();
        lsVTKWriter(mesh, "sphere-" + std::to_string(counter) + ".vtk").apply();

        //lsToMesh<double, D>(newLayer, mesh).apply();
        //lsVTKWriter(mesh, "LS-" + std::to_string(counter) + ".vtk").apply();

        ++counter;
    }
*/
    return 0;


}