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



#include <lsEikonalExpand.hpp>



#include <../dev/derivatives.hpp>

#include <omp.h>

//#include <hrleSparseStarIterator.hpp>

//#include <hrleCartesianPlaneIterator.hpp>


//____________testing not necessary_________________

#include <chrono>

#include <../dev/lsEikonalExpandTest.hpp>

//____________testing end___________________________

constexpr int D = 3 ;
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

lsSmartPointer<lsDomain<double, D>> makePlane(double gridDelta, std::vector<NumericType>& planeNormal){

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

lsSmartPointer<lsDomain<double, D>> makeTrench(double gridDelta, std::vector<NumericType>& planeNormal){

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

void create_output(lsSmartPointer<lsDomain<double, D>> levelSet,
                   std::string output_name){

    std::vector<NumericType> meanCurvatureGeneralFormula;

    std::vector<NumericType> meanCurvatureGeneralFormulaBig;

    std::vector<NumericType> meanVariationOfNormal;

    std::vector<NumericType> meanCurvatureGeneralFormulaBigBias;

    std::vector<NumericType> meanCurvatureSD1;

    std::vector<NumericType> meanCurvature2;

    std::vector<NumericType> meanCurvatureSD2;

    std::vector<NumericType> meanCurveShapeBias;


    NumericType gridDelta = levelSet->getGrid().getGridDelta();

    auto grid = levelSet->getGrid();

    typename lsDomain<NumericType, D>::DomainType &domain = levelSet->getDomain();

    double pointsPerSegment =
    double(2 * levelSet->getDomain().getNumberOfPoints()) /
    double(levelSet->getLevelSetWidth());

    std::vector<std::vector<NumericType>> meanCurvatureReserve(levelSet->getNumberOfSegments());


#pragma omp parallel num_threads((levelSet)->getNumberOfSegments())
        {
            int p = 0;
#ifdef _OPENMP
            p = omp_get_thread_num();
#endif
            curvaturShapeDerivatives1<NumericType, D> shape1(gridDelta);

            std::vector<NumericType> &meanCurveSegment = meanCurvatureReserve[p];
            meanCurveSegment.reserve(pointsPerSegment);
           
            hrleVectorType<hrleIndexType, D> startVector =
            (p == 0) ? grid.getMinGridPoint()
                : domain.getSegmentation()[p - 1];

            hrleVectorType<hrleIndexType, D> endVector =
            (p != static_cast<int>(domain.getNumberOfSegments() - 1))
                ? domain.getSegmentation()[p]
                : grid.incrementIndices(grid.getMaxGridPoint());


            for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
                    neighborIt(levelSet->getDomain(), startVector);
                    neighborIt.getIndices() < endVector; neighborIt.next()) {
                   
                //TODO: add value for EUCLID normalization!!!
                if (!neighborIt.getCenter().isDefined() || std::abs(neighborIt.getCenter().getValue()) > 0.5) {
                       continue;
                } 

                meanCurveSegment.push_back(shape1(neighborIt));
           
            }

        }

        meanCurvatureSD1.reserve(levelSet->getNumberOfPoints()); 

        for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i){ 
            meanCurvatureSD1.insert(meanCurvatureSD1.end(), meanCurvatureReserve[i].begin(),
             meanCurvatureReserve[i].end());
        }

        std::vector<std::vector<NumericType>> generalFormulaReserve(levelSet->getNumberOfSegments());

        std::vector<std::vector<NumericType>> variationOfGradReserve(levelSet->getNumberOfSegments());

        std::vector<std::vector<NumericType>> generalFormulaBigReserve(levelSet->getNumberOfSegments());
        
        std::vector<std::vector<NumericType>> shape2Reserve(levelSet->getNumberOfSegments());

        std::vector<std::vector<NumericType>> generalFormulaBigBiasReserve(levelSet->getNumberOfSegments());

        std::vector<std::vector<NumericType>> shapeBiasReserve(levelSet->getNumberOfSegments());


#pragma omp parallel num_threads((levelSet)->getNumberOfSegments())
        {
            int p = 0;
#ifdef _OPENMP
            p = omp_get_thread_num();
#endif

            curvaturGeneralFormula<NumericType, D> generalFormula(gridDelta);

            variationOfNormals<NumericType, D> variationOfGrad(gridDelta);

            curvaturGeneralFormulaBigStencil <NumericType, D> generalFormulaBig(gridDelta);

            curvaturShapeDerivatives2<NumericType, D> shape2(gridDelta);

            curvaturGeneralFormulaBigStencilBias <NumericType, D> generalFormulaBigBias(gridDelta);

            curvaturShapeBias<NumericType, D> shapeBias(gridDelta);

            std::vector<NumericType> &generalFormulaSegment = generalFormulaReserve[p];
            generalFormulaSegment.reserve(pointsPerSegment);

            std::vector<NumericType> &variationOfGradSegment = variationOfGradReserve[p];
            variationOfGradSegment.reserve(pointsPerSegment);

            std::vector<NumericType> &generalFormulaBigSegment = generalFormulaBigReserve[p];
            generalFormulaBigSegment.reserve(pointsPerSegment);

            std::vector<NumericType> &shape2Segment = shape2Reserve[p];
            shape2Segment.reserve(pointsPerSegment);

            std::vector<NumericType> &generalFormulaBigBiasSegment = generalFormulaBigBiasReserve[p];
            generalFormulaBigBiasSegment.reserve(pointsPerSegment);

            std::vector<NumericType> &shapeBiasSegment = shapeBiasReserve[p];
            shapeBiasSegment.reserve(pointsPerSegment);

            hrleVectorType<hrleIndexType, D> startVector =
            (p == 0) ? grid.getMinGridPoint()
                : domain.getSegmentation()[p - 1];

            hrleVectorType<hrleIndexType, D> endVector =
            (p != static_cast<int>(domain.getNumberOfSegments() - 1))
                ? domain.getSegmentation()[p]
                : grid.incrementIndices(grid.getMaxGridPoint());


            for (hrleCartesianPlaneIterator<typename lsDomain<NumericType, D>::DomainType>
                 neighborIt(levelSet->getDomain(), startVector, 1);
                 neighborIt.getIndices() < endVector; neighborIt.next()) {
                   
                //TODO: add value for Euler normalization!!!
                if (!neighborIt.getCenter().isDefined() || std::abs(neighborIt.getCenter().getValue()) > 0.5) {
                       continue;
                } 

                generalFormulaSegment.push_back(generalFormula(neighborIt));

                variationOfGradSegment.push_back(variationOfGrad(neighborIt));

                generalFormulaBigSegment.push_back(generalFormulaBig(neighborIt));

                shape2Segment.push_back(shape2(neighborIt));

                generalFormulaBigBiasSegment.push_back(generalFormulaBigBias(neighborIt));

                shapeBiasSegment.push_back(shapeBias(neighborIt));
           
            }

        }

        for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i){ 
            meanCurvatureGeneralFormula.insert(meanCurvatureGeneralFormula.end(), generalFormulaReserve[i].begin(),
             generalFormulaReserve[i].end());
        }

        for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i){ 
            meanVariationOfNormal.insert(meanVariationOfNormal.end(), variationOfGradReserve[i].begin(),
             variationOfGradReserve[i].end());
        }

        for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i){ 
            meanCurvatureGeneralFormulaBig.insert(meanCurvatureGeneralFormulaBig.end(), generalFormulaBigReserve[i].begin(),
             generalFormulaBigReserve[i].end());
        }

        for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i){ 
            meanCurvatureSD2.insert(meanCurvatureSD2.end(), shape2Reserve[i].begin(),
             shape2Reserve[i].end());
        }

        for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i){ 
            meanCurvatureGeneralFormulaBigBias.insert(meanCurvatureGeneralFormulaBigBias.end(), generalFormulaBigBiasReserve[i].begin(),
             generalFormulaBigBiasReserve[i].end());
        }

        for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i){ 
            meanCurveShapeBias.insert(meanCurveShapeBias.end(), shapeBiasReserve[i].begin(),
             shapeBiasReserve[i].end());
        }
    //TODO: check if gd has to be 
    auto narrowband = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;

    if(levelSet->getLevelSetNormalization() == lsNormalizations::MANHATTEN){
      lsToMesh<NumericType, D>(levelSet, narrowband, true, true).apply();
    }else{
      lsToMesh<NumericType, D>(levelSet, narrowband, true, true, gridDelta).apply();
    }
    

    //narrowband.insertNextVectorData(normal, "Normal");
    narrowband->insertNextScalarData(meanCurvatureSD1, "shape operator small stencil");
    narrowband->insertNextScalarData(meanCurvatureGeneralFormula, "general formula");
    narrowband->insertNextScalarData(meanVariationOfNormal, "variatoin of normal");
    narrowband->insertNextScalarData(meanCurvatureGeneralFormulaBig, "general formula big stencil");
    narrowband->insertNextScalarData(meanCurvatureSD2, "shape operator big stencil");
    narrowband->insertNextScalarData(meanCurvatureGeneralFormulaBigBias, "general formula big stencil bias");
    narrowband->insertNextScalarData(meanCurveShapeBias, "shape operator Bias");

  
    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , output_name ).apply();


  }





int main(int argc, char* argv[]) {

    lsNormalizations normalization = lsNormalizations::EUCLID;

    omp_set_num_threads(4);

    NumericType gridDelta = 0.5;

    NumericType radius = 10.;

    std::string outputName = "curveOutput";



    if(argc == 5){
        if(std::string(argv[1]) != "euclid"){
           normalization = lsNormalizations::MANHATTEN;
        }

        gridDelta = std::stod(argv[2]);

        radius = std::stod(argv[3]);

        outputName = std::string(argv[4]);
        std::cout << "#SOMEchanges" << std::endl;
    }else{
      std::cout << "#nochanges" << std::endl;
    }


    std::vector<lsSmartPointer<lsDomain<double, D>>> levelSets;

    lsSmartPointer<lsDomain<double, D>> levelSet = makeSphere(gridDelta, radius, normalization);

    //std::vector<NumericType> normal = {0., 0., 1.};

    //lsSmartPointer<lsDomain<double, D>> levelSet = makeTrench(gridDelta, normal);

 
    levelSets.push_back(levelSet);  


    std::cout << "Fast Marching..." << std::endl;
    if(normalization == lsNormalizations::EUCLID){
      lsEikonalExpandTest<NumericType, D> expander(levelSets.back(), 5);
      expander.apply(); 
    }else{
      lsExpand<NumericType, D>(levelSets.back(), 7).apply();
    }
   
    std::cout << "Calculating Curvatures..." << std::endl;

    create_output(levelSets.back(), outputName);

    std::cout << "Finished" << std::endl;

    return 0;
}





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




