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



//____________testing not necessary_________________

#include <chrono>

#include <lsCalculateNormalVectors.hpp>

#include <../dev/lsEikonalExpandTest.hpp>

//____________testing end___________________________

constexpr int D = 3 ;
typedef double NumericType;

typedef typename lsDomain<NumericType, D>::DomainType hrleDomainType;

void timingTests(lsSmartPointer<lsDomain<double, D>> levelSet,
  std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & lsPoints){

    hrleSparseBoxIterator<hrleDomain<NumericType, D>> neighborIterator(levelSet->getDomain(), 1);

    hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType> neighborStarIterator(levelSet->getDomain());
       
    NumericType gridDelta = 1.;//levelSet.getGrid().getGridDelta();

    variationOfNormals<NumericType, D> variationOfGrad(gridDelta);

    curvaturShapeDerivatives1<NumericType, D> shape1(gridDelta);

    curvaturShapeDerivatives2<NumericType, D> shape2(gridDelta);

    curvaturGeneralFormula<NumericType, D> generalFormula(gridDelta);

    curvaturShapeBias<NumericType, D> shapeBias(gridDelta);

    curvaturGeneralFormulaBigStencil <NumericType, D> generalFormulaBig(gridDelta);

    curvaturGeneralFormulaBigStencilBias <NumericType, D> generalFormulaBigBias(gridDelta);

    std::cout << "Doing Timing Test..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now(); 

    auto stop = std::chrono::high_resolution_clock::now(); 

    int runs = 50;

    double time = 0.;
    std::cout << "meanCurvatureGeneralFormula" <<std::endl;
    for(int i = 0;i < runs; i ++){
      start = std::chrono::high_resolution_clock::now(); 
    for(hrleConstSparseIterator<hrleDomainType> centerIt(levelSet->getDomain());
      !centerIt.isFinished(); ++centerIt){
   // for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
    //      neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
     //     neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

        //auto &centerIt = neighborIt.getCenter();

              if (!centerIt.isDefined() || (lsPoints.find(centerIt.getStartIndices()) == lsPoints.end())) {
                continue;
              } 

              // move neighborIterator to current position
              neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

              generalFormula(neighborIterator);             

      }
      stop = std::chrono::high_resolution_clock::now(); 

      time += std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

      std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << ", "; 
    }

    std::cout << std::endl << (time/((double)runs)) << std::endl;
    time = 0.;

    std::cout << "meanCurvatureGeneralFormulaBigStencil" <<std::endl;
    for(int i = 0;i < runs; i ++){
      start = std::chrono::high_resolution_clock::now(); 
      for(hrleConstSparseIterator<hrleDomainType> centerIt(levelSet->getDomain());
        !centerIt.isFinished(); ++centerIt){
      // for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
      //      neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
      //     neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

              //auto &centerIt = neighborIt.getCenter();
              if (!centerIt.isDefined() || (lsPoints.find(centerIt.getStartIndices()) == lsPoints.end())) {
                continue;
              } 

              // move neighborIterator to current position
              neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

              generalFormulaBig(neighborIterator);             

      }
      stop = std::chrono::high_resolution_clock::now(); 

      time += std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

      std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << ", "; 
    }

    std::cout << std::endl << (time/((double)runs)) << std::endl;
    time = 0.;

    std::cout << "meanCurvatureGeneralFormulaBigBias" <<std::endl;
    for(int i = 0;i < runs; i ++){
      start = std::chrono::high_resolution_clock::now(); 
      for(hrleConstSparseIterator<hrleDomainType> centerIt(levelSet->getDomain());
        !centerIt.isFinished(); ++centerIt){
      // for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
      //      neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
      //     neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

        //auto &centerIt = neighborIt.getCenter();
              if (!centerIt.isDefined() || (lsPoints.find(centerIt.getStartIndices()) == lsPoints.end())) {
                continue;
              } 

              // move neighborIterator to current position
              neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

              generalFormulaBigBias(neighborIterator);             

      }
      stop = std::chrono::high_resolution_clock::now(); 

      time += std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

      std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << ", "; 
    }

    std::cout << std::endl << (time/((double)runs)) << std::endl;
    time = 0.;


    std::cout << "VariationOfGradient" <<std::endl;
    for(int i = 0;i < runs; i ++){
      start = std::chrono::high_resolution_clock::now(); 
      for(hrleConstSparseIterator<hrleDomainType> centerIt(levelSet->getDomain());
        !centerIt.isFinished(); ++centerIt){
        // for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
        //      neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
        //     neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

        //auto &centerIt = neighborIt.getCenter();
              if (!centerIt.isDefined() || (lsPoints.find(centerIt.getStartIndices()) == lsPoints.end())) {
                continue;
              } 

              // move neighborIterator to current position
              neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

              variationOfGrad(neighborIterator);             

      }
      stop = std::chrono::high_resolution_clock::now(); 

      time += std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

      std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << ", "; 
    }

    std::cout << std::endl << (time/((double)runs)) << std::endl;
    time = 0.;

    std::cout << "ShapeBias" <<std::endl;

    for(int i = 0;i < runs; i ++){
      start = std::chrono::high_resolution_clock::now(); 
      for(hrleConstSparseIterator<hrleDomainType> centerIt(levelSet->getDomain());
        !centerIt.isFinished(); ++centerIt){
    // for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
      //      neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
      //     neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

        //auto &centerIt = neighborIt.getCenter();
              if (!centerIt.isDefined() || (lsPoints.find(centerIt.getStartIndices()) == lsPoints.end())) {
                continue;
              } 

              // move neighborIterator to current position
              neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

              shapeBias(neighborIterator);             

      }
      stop = std::chrono::high_resolution_clock::now(); 

      time += std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

      std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << ", "; 
    }

    std::cout << std::endl << (time/((double)runs)) << std::endl;
    time = 0.;

    std::cout << "Shape1" <<std::endl;

    for(int i = 0;i < runs; i ++){
      start = std::chrono::high_resolution_clock::now(); 
      for(hrleConstSparseIterator<hrleDomainType> centerIt(levelSet->getDomain());
        !centerIt.isFinished(); ++centerIt){
        // for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
        //      neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
        //     neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

          //auto &centerIt = neighborIt.getCenter();
              if (!centerIt.isDefined() || (lsPoints.find(centerIt.getStartIndices()) == lsPoints.end())) {
                continue;
              } 

              // move neighborStarIterator to current position
              neighborStarIterator.goToIndicesSequential(centerIt.getStartIndices());

              shape1(neighborStarIterator);             

      }
      stop = std::chrono::high_resolution_clock::now(); 

      time += std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

      std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << ", "; 
    }

    std::cout << std::endl << (time/((double)runs)) << std::endl;
    time = 0.;

    std::cout << "Shape2" <<std::endl;

    for(int i = 0;i < runs; i ++){
      start = std::chrono::high_resolution_clock::now(); 
      for(hrleConstSparseIterator<hrleDomainType> centerIt(levelSet->getDomain());
        !centerIt.isFinished(); ++centerIt){
      // for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
      //      neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
      //     neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

        //auto &centerIt = neighborIt.getCenter();
             
              if (!centerIt.isDefined() || (lsPoints.find(centerIt.getStartIndices()) == lsPoints.end())) {
                continue;
              } 

              // move neighborIterator to current position
              neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

              shape2(neighborIterator);             

      }
      stop = std::chrono::high_resolution_clock::now(); 

      time += std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

      std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << ", "; 
    }

    std::cout << std::endl << (time/((double)runs)) << std::endl;
    time = 0.; 

}



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

lsDomain<double, D> makePlane(double gridDelta, std::vector<NumericType>& planeNormal){

      std::cout << "creating Plane..." << std::endl;

    double extent = 5;
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

    auto levelSet =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);


    std::vector<NumericType> origin = {0., 0., 0.};


    {
      auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
      lsMakeGeometry<double, D>(levelSet, plane).apply();
    }

    return levelSet;


}

lsDomain<double, D> makeTrench(double gridDelta, std::vector<NumericType>& planeNormal){

    std::cout << "creating trench..." << std::endl;

    double extent = 50;
    double bounds[2 * D] = {-extent, extent, -extent, extent};
    if (D == 3) {
        bounds[4] = -extent;
        bounds[5] = extent;
    }

    typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
    for (unsigned i = 0; i < D - 1; ++i) {
        boundaryCons[i] =
            lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    }
    boundaryCons[D - 1] =
        lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

    auto levelSet =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

    std::vector<NumericType> origin = {0., 0., 0.};

    {
      auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
      lsMakeGeometry<double, D>(levelSet, plane).apply();
    }

    {
        // create layer used for booling
        std::cout << "Creating box..." << std::endl;

        auto trench = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

        if(D == 3){
          double minCorner[3] = {-extent - 1, -extent / 4., -15.};
          double maxCorner[3] = {extent + 1, extent / 4., 1.0};
          auto box = lsSmartPointer<lsBox<double, D>>::New(minCorner, maxCorner);
          lsMakeGeometry<double, D>(trench, box).apply();
        }else{
          double minCorner[2] = {-extent / 4., -15.};
          double maxCorner[2] = {extent / 4., 1.0};
          auto box = lsSmartPointer<lsBox<double, D>>::New(minCorner, maxCorner);
          lsMakeGeometry<double, D>(trench, box).apply();
        }
        //lsMakeGeometry<double, D>(trench, lsBox<double, D>(minCorner, maxCorner))
        //    .apply();


        // Create trench geometry
        std::cout << "Booling trench..." << std::endl;
        lsBooleanOperation<double, D>(levelSet, trench,
                                      lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
    }

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
  std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & lsPoints,
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
        if (!centerIt.isDefined() || std::abs(centerIt.getValue()) > gridDelta) {
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
    lsToMesh<NumericType, D>(levelSet, narrowband, true, true).apply(lsPoints);
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





int main() {


    omp_set_num_threads(1);

    NumericType gridDelta = 0.5;

    //______________________________First____________________________________________________________________

    //std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints;

    std::vector<lsSmartPointer<lsDomain<double, D>>> levelSets;

      //90 grad
    //std::vector<NumericType> planeNormal = {0. , 0. , 1.};

    //22.5 grad
    //planeNormal = {0.0, 0.38268343, 0.92387953};
    // additionalParams.push_back(planeNormal);

    //45 grad
    //planeNormal = {0.0, 0.70710678, 0.70710678};
    //additionalParams.push_back(planeNormal);

    //67.5 grad
    //planeNormal ={0.0, 0.92387953, 0.38268343};
    //additionalParams.push_back(planeNormal);

    //lsDomain<NumericType,D> levelSet = makePlane(gridDelta, planeNormal);

    NumericType radius = 100.;

    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> activePoints;

    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> narrowPoints;

    lsSmartPointer<lsDomain<double, D>> levelSet = makeSphere(gridDelta, radius, activePoints, narrowPoints);

    //lsDomain<NumericType,D> levelSet = makeTrench(gridDelta, planeNormal);





    levelSets.push_back(levelSet);  

/*
    int order = 1;
    lsExpand<NumericType, D>(levelSets.back(), 2 * (order + 2) + 1).apply();

    std::cout << "Calculating Curvatures..." << std::endl;

    create_output_Manhatten(levelSets.back(),  "final_output_Manhatten");
*/


/*
    std::cout << "Converting..." << std::endl;

    lsConvertEuclid<NumericType, D>  converter(levelSets.back());

    converter.apply();



    //get the active grid points of the level set
    activePoints = converter.getActivePoints();

    std::cout << activePoints.size() << std::endl;

*/

    auto start = std::chrono::high_resolution_clock::now(); 

    std::cout << "Fast Marching..." << std::endl;

    //lsEikonalExpand<NumericType, D> expander(levelSets.back(), narrowPoints);

    lsEikonalExpandTest<NumericType, D> expander(levelSets.back());

    //lsExpandSphere<NumericType, D> expander(levelSets.back(), activePoints, radius);

    expander.apply(); 

    auto stop = std::chrono::high_resolution_clock::now(); 

    std::cout << "time for FMM: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << std::endl; 

    std::cout << "Calculating Curvatures..." << std::endl;

    create_output_Euklid(levelSets.back(), narrowPoints, "/media/sf_shared/narrowband");


    //timingTests(levelSets.back(), activePoints);


    std::cout << "Finished" << std::endl;

    return 0;
}




