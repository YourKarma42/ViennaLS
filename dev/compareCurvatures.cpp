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
#include <lsEiknoalExpand.hpp>

//#include <hrleSparseBoxIterator.hpp>
//#include <hrleVectorType.hpp>



#include <lsCalculateCurvatures.hpp>

#include <omp.h>

#include <../dev/derivatives.hpp>



//____________testing not necessary_________________

#include <chrono>

#include <lsCalculateNormalVectors.hpp>

//____________testing end___________________________

constexpr int D = 3;
typedef double NumericType;



lsDomain<double, D> makeSphere(double gridDelta, double radius){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0., 0., 0.};
    
    lsDomain<double,D> levelSet(gridDelta);

    lsMakeGeometry<double, D>(levelSet, lsSphere<double, D>(origin, radius)).apply();


    return levelSet;

}

void create_output(lsDomain<NumericType,D> & levelSet,
  std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & lsPoints,
  std::string output_name){

    hrleSparseBoxIterator<hrleDomain<NumericType, D>> neighborIterator(levelSet.getDomain(), 2);

    std::vector<std::array<NumericType, 3>> normal;

    std::vector<std::array<NumericType, 3>> secondOrderDerivatives1;

    std::vector<std::array<NumericType, 3>> secondOrderDerivatives2;



    std::vector<NumericType> grad1;

    std::vector<NumericType> meanCurvatureGeneralFormula;

    std::vector<NumericType> meanCurvatureSD1;

    std::vector<NumericType> grad2;

    std::vector<NumericType> meanCurvature2;

    std::vector<NumericType> meanCurvatureSD2;

    std::vector<NumericType> meanCurvature3;

    std::vector<NumericType> meanCurvatureNew;

    std::vector<NumericType> curveTest;

    std::vector<NumericType> meanCurveShapeBias;


    //std::vector<std::array<NumericType, D>> normal;


    NumericType gridDelta = 1.;//levelSet.getGrid().getGridDelta();

    variationOfNormals<NumericType, D> variationOfGrad(gridDelta);

    curvaturShapeDerivatives1<NumericType, D> shape1(gridDelta);

    curvaturShapeDerivatives2<NumericType, D> shape2(gridDelta);

    curvaturGeneralFormula<NumericType, D> generalFormula(gridDelta);

    curvaturShapeBias<NumericType, D> shapeBias(gridDelta);

    curvaturTest<NumericType, D> test(gridDelta);




   
   //     ? passedlsDomain.getDomain().getSegmentation()[p]
   //     : grid.incrementIndices(grid.getMaxGridPoint());

    for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
          neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
          neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

        auto &centerIt = neighborIt.getCenter();
        if (!centerIt.isDefined() || (lsPoints.find(centerIt.getStartIndices()) == lsPoints.end())) {
          continue;
        } 

        // move neighborIterator to current position
        neighborIterator.goToIndicesSequential(centerIt.getStartIndices());


        meanCurvatureNew.push_back(variationOfGrad(neighborIterator));

        meanCurvatureSD1.push_back(shape1(neighborIterator));

        meanCurvatureSD2.push_back(shape2(neighborIterator));

        meanCurvatureGeneralFormula.push_back(generalFormula(neighborIterator));

        meanCurveShapeBias.push_back(shapeBias(neighborIterator));

        curveTest.push_back(test(neighborIterator));


    }

    lsMesh narrowband;
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(levelSet, narrowband, true, true).apply(lsPoints);
    //lsPoints


    //narrowband.insertNextVectorData(normal, "Normal");
    narrowband.insertNextScalarData(meanCurvatureSD1, "shape operator small stencil");
    narrowband.insertNextScalarData(meanCurvatureSD2, "shape operator big stencil");
    narrowband.insertNextScalarData(meanCurvatureGeneralFormula, "general formula");
    narrowband.insertNextScalarData(curveTest, "Curvature Test");
    narrowband.insertNextScalarData(meanCurvatureNew, "mean new");
    narrowband.insertNextScalarData(meanCurveShapeBias, "shape operator Bias");

  
    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , output_name ).apply();

    std::cout << "Doing Timing Test..." << std::endl;

    //________________________________________Development area_______________________________________________________________________________________________________________________    

   /* auto start = std::chrono::high_resolution_clock::now(); 

    auto stop = std::chrono::high_resolution_clock::now(); 

    int runs = 50;

    double time = 0.;
    std::cout << "meanCurvatureGeneralFormula" <<std::endl;
    for(int i = 0;i < runs; i ++){
      start = std::chrono::high_resolution_clock::now(); 
      for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
            neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
            neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

              auto &centerIt = neighborIt.getCenter();
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

    std::cout << "VariationOfGradient" <<std::endl;
    for(int i = 0;i < runs; i ++){
      start = std::chrono::high_resolution_clock::now(); 
      for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
            neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
            neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

              auto &centerIt = neighborIt.getCenter();
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
      for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
            neighborIt1(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
            neighborIt1.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt1.next()) {

              auto &centerIt = neighborIt1.getCenter();
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
      for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
            neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
            neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

              auto &centerIt = neighborIt.getCenter();
              if (!centerIt.isDefined() || (lsPoints.find(centerIt.getStartIndices()) == lsPoints.end())) {
                continue;
              } 

              // move neighborIterator to current position
              neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

              shape1(neighborIterator);             

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
      for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
            neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
            neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

              auto &centerIt = neighborIt.getCenter();
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
    time = 0.; */


}




int main() {


    omp_set_num_threads(4);

    NumericType gridDelta = 0.25;

    //______________________________First____________________________________________________________________

    //std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints;

    std::vector<lsDomain<NumericType, D> *> levelSets;

    lsDomain<NumericType,D> levelSet = makeSphere(gridDelta, 10.);

    levelSets.push_back(&levelSet);  

    std::cout << "Converting..." << std::endl;

    lsConvertEuclid<NumericType, D>  converter(*(levelSets.back()));

    converter.apply();

    //get the active grid points of the level set
    auto activePoints = converter.getActivePoints();

    std::cout << "Fast Marching..." << std::endl;

    lsEikonalExpand<NumericType, D> expander(*(levelSets.back()), activePoints);

    expander.apply(); 

    create_output(*(levelSets.back()), activePoints, "final_output");


    std::cout << "Finished" << std::endl;

    //______________________________Second____________________________________________________________________
 
   /* std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints1;

    std::vector<lsDomain<NumericType, D> *> levelSets1;

    lsDomain<NumericType,D> levelSet1 = makeSphere(gridDelta, 15.);

    lsMesh diskmesh1;

    lsToDiskMesh<NumericType, D>(levelSet1, diskmesh1).apply();

    lsVTKWriter(diskmesh1, lsFileFormatEnum::VTU , "SphereDiskMesh" ).apply();

    levelSets1.push_back(&levelSet1);  

    lsDomain<NumericType,D> newLevelSet1 =  ConvertLS(*(levelSets1.back()), lsPoints1);

    testFMM(& newLevelSet1, lsPoints1);

    create_output(newLevelSet1, lsPoints1, "radius_50");*/


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

