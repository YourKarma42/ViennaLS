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

//#include <lsSmartPointer.hpp>


#include <unordered_set>

#include <filesystem>



#include <lsConvertEuclid.hpp>
#include <lsEikonalExpand.hpp>

#include "hrleTestIterator.hpp"


#include <omp.h>


#include <../dev/lsCurvatureCalculator.hpp>


//____________testing not necessary_________________

#include <chrono>

#include <lsCalculateNormalVectors.hpp>

//____________testing end___________________________

constexpr int D = 3;
typedef double NumericType;

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


int main(int argc, char* argv[]) {

    std::stringstream csvOutput;

    int numThreads = 1;

    int numberOfRuns = 10;

    std::vector<double> timings;

    if(argc != 1){
        numThreads = std::stoi(argv[1]);
    }

    std::cout << "running program with " << numThreads << " threads" << std::endl;

    csvOutput << numberOfRuns << std::endl;

    csvOutput << numThreads << std::endl;


    omp_set_num_threads(numThreads);

    NumericType gridDelta = 0.25;

    std::cout << "grid delta is " << gridDelta << std::endl; 

    auto start = std::chrono::high_resolution_clock::now(); 

    auto stop = std::chrono::high_resolution_clock::now(); 

    std::vector<lsSmartPointer<lsDomain<double, D>>> levelSets;

    NumericType radius = 100.;

    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> narrowPoints;

    lsSmartPointer<lsDomain<double, D>> levelSet = makeSphere(gridDelta, radius, narrowPoints);

    //lsDomain<NumericType,D> levelSet = makeTrench(gridDelta);

    levelSets.push_back(levelSet);  

    std::cout << "geometry is Sphere with radius: " << radius << std::endl; 

    std::cout << "FMM..." << std::endl;

    lsEikonalExpand<NumericType, D> expander(levelSets.back(), narrowPoints);

    expander.apply(); 


//_______________________Prepare stuff___________________________
    std::vector<lsSmartPointer<lsInternal::baseCurvature<NumericType,D>>> vecDerivatives;
    //std::vector<baseDerivative<NumericType,D>*> vecDerivatives;
    std::vector<std::string> names;

    vecDerivatives.push_back(lsSmartPointer<lsInternal::curvaturGeneralFormula<NumericType,D>>::New(gridDelta));
    names.push_back("General Formula");

    vecDerivatives.push_back(lsSmartPointer<lsInternal::curvaturGeneralFormulaBigStencil<NumericType,D>>::New(gridDelta));
    names.push_back("General Formula Big Stencil");

    vecDerivatives.push_back(lsSmartPointer<lsInternal::variationOfNormals<NumericType,D>>::New(gridDelta));
    names.push_back("Variation of Normals");

    vecDerivatives.push_back(lsSmartPointer<lsInternal::curvaturGeneralFormulaBigStencilBias<NumericType,D>>::New(gridDelta));
    names.push_back("General Formula Bias");

    vecDerivatives.push_back(lsSmartPointer<lsInternal::curvaturShapeDerivatives2<NumericType,D>>::New(gridDelta));
    names.push_back("Shape Operator Derivatives 2");

    vecDerivatives.push_back(lsSmartPointer<lsInternal::curvaturShapeBias<NumericType,D>>::New(gridDelta));
    names.push_back("Shape Operator Bias");





    NumericType tmp=0.;

//______________________________________________________Start______________________________________________________________

    auto namesItr = names.begin();

    for(auto curve: vecDerivatives){

        std::cout << *namesItr << std::endl;  

        csvOutput << *namesItr << std::endl;

        namesItr++;
   
        for(int i = 0; i < numberOfRuns; i++){

            lsInternal::curvaturGeneralFormula<NumericType, D> generalFormula(gridDelta);

            auto grid = levelSets.back()->getGrid();

            typename lsDomain<NumericType, D>::DomainType &domain = levelSets.back()->getDomain();

            double pointsPerSegment =
            double(2 * levelSets.back()->getDomain().getNumberOfPoints()) /
            double(levelSets.back()->getLevelSetWidth());

            std::vector<std::vector<NumericType>> meanCurvatureReserve(levelSets.back()->getNumberOfSegments());


            start = std::chrono::high_resolution_clock::now(); 


#pragma omp parallel num_threads((levelSets.back())->getNumberOfSegments())
        {
            int p = 0;
#ifdef _OPENMP
            p = omp_get_thread_num();
#endif

            std::vector<NumericType> &meanCurveSegment = meanCurvatureReserve[p];
            meanCurveSegment.reserve(pointsPerSegment);
            hrleCartesianPlaneIterator<typename lsDomain<NumericType, D>::DomainType> neighborIt(domain, 1);

            hrleVectorType<hrleIndexType, D> startVector =
            (p == 0) ? grid.getMinGridPoint()
                : domain.getSegmentation()[p - 1];

            hrleVectorType<hrleIndexType, D> endVector =
            (p != static_cast<int>(domain.getNumberOfSegments() - 1))
                ? domain.getSegmentation()[p]
                : grid.incrementIndices(grid.getMaxGridPoint());


            for(hrleSparseIterator<typename lsDomain<NumericType, D>::DomainType> it(
                domain, startVector);
                it.getStartIndices() < endVector; ++it){

                if (!it.isDefined() || std::abs(it.getValue()) > gridDelta) {
                    continue;
                }
                
                neighborIt.goToIndices(it.getStartIndices());
                //curve->calcDerivatives(neighborIt);
                meanCurveSegment.push_back(curve->getMeanCurvature(neighborIt));
                //tmp = curve->getGaussianCurvature();          
            }

        }

        std::vector<NumericType> meanCurvature;

        meanCurvature.reserve(levelSets.back()->getNumberOfPoints()); 

        for (unsigned i = 0; i < levelSets.back()->getNumberOfSegments(); ++i){ 
            meanCurvature.insert(meanCurvature.end(), meanCurvatureReserve[i].begin(),
             meanCurvatureReserve[i].end());
        }

        stop = std::chrono::high_resolution_clock::now();

        csvOutput << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << ";";

        }
        csvOutput << std::endl;
    }



//______________________________________________________End______________________________________________________________

//______________________________________________________Start______________________________________________________________

    std::cout << "Shape Operator Small Stencil" << std::endl;  

    csvOutput << "Shape Operator Small Stencil" << std::endl;

    for(int i = 0; i < numberOfRuns; i++){

        lsInternal::curvaturShapeDerivatives1<NumericType, D> shapeOperator(gridDelta);

        auto grid = levelSets.back()->getGrid();

        typename lsDomain<NumericType, D>::DomainType &domain = levelSets.back()->getDomain();

        double pointsPerSegment =
        double(2 * levelSets.back()->getDomain().getNumberOfPoints()) /
        double(levelSets.back()->getLevelSetWidth());

        std::vector<std::vector<NumericType>> meanCurvatureReserve(levelSets.back()->getNumberOfSegments());


        start = std::chrono::high_resolution_clock::now(); 


#pragma omp parallel num_threads((levelSets.back())->getNumberOfSegments())
        {
            int p = 0;
#ifdef _OPENMP
            p = omp_get_thread_num();
#endif

            std::vector<NumericType> &meanCurveSegment = meanCurvatureReserve[p];
            meanCurveSegment.reserve(pointsPerSegment);
            hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType> neighborIt(domain);

            hrleVectorType<hrleIndexType, D> startVector =
            (p == 0) ? grid.getMinGridPoint()
                : domain.getSegmentation()[p - 1];

            hrleVectorType<hrleIndexType, D> endVector =
            (p != static_cast<int>(domain.getNumberOfSegments() - 1))
                ? domain.getSegmentation()[p]
                : grid.incrementIndices(grid.getMaxGridPoint());


            for(hrleSparseIterator<typename lsDomain<NumericType, D>::DomainType> it(
                domain, startVector);
                it.getStartIndices() < endVector; ++it){

                if (!it.isDefined() || std::abs(it.getValue()) > gridDelta) {
                    continue;
                }
                
                neighborIt.goToIndices(it.getStartIndices());
                //std::cout << "TÜ" << std::endl;
                //shapeOperator.calcDerivatives(neighborIt);
                //meanCurveSegment.push_back(shapeOperator.getGaussianCurvature());
                meanCurveSegment.push_back(shapeOperator.getMeanCurvature(neighborIt));

            
            }

        }

        std::vector<NumericType> meanCurvature;

        meanCurvature.reserve(levelSets.back()->getNumberOfPoints()); 

        for (unsigned i = 0; i < levelSets.back()->getNumberOfSegments(); ++i){ 
            meanCurvature.insert(meanCurvature.end(), meanCurvatureReserve[i].begin(),
             meanCurvatureReserve[i].end());
        }

        stop = std::chrono::high_resolution_clock::now();

        csvOutput << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << ";";

    }

    csvOutput << std::endl;
    std::cout << "blub" << std::endl;

//______________________________________________________End______________________________________________________________

    //std::filesystem::create_directory("timings");

    std::ofstream output;

    output.open("timings/timingsCurvature" + std::to_string(numThreads) + ".csv");

    output << csvOutput.rdbuf();

    output.close();

    //std::cout << csvOutput.str() << std::endl;

    std::cout << "Finished" << std::endl;

    return 0;
}




