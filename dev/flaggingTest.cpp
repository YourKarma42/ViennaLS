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

#include <lsReader.hpp>

#include <lsConvertEuclid.hpp>
#include <../dev/lsEikonalExpandTest.hpp>

#include <../dev/silvacoLikeFlagging.hpp>
#include <../dev/myFlagging.hpp>


#include <omp.h>



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
          double minCorner[3] = {-extent / 4., -extent - 1, -50.};
          double maxCorner[3] = { extent / 4., extent + 1, 1.0};
          auto box = lsSmartPointer<lsBox<double, D>>::New(minCorner, maxCorner);
          lsMakeGeometry<double, D>(trench, box).apply();
        }else{
          double minCorner[2] = {-extent / 4., -50.};
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

    levelSet->getDomain().segment();

    return levelSet;
    
}

lsSmartPointer<lsDomain<double, D>> makeSphere(double gridDelta, double radius){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0.0, 0.0, 0.0};
    
    auto levelSet =
        lsSmartPointer<lsDomain<double, D>>::New(gridDelta);

    lsMakeGeometry<double, D>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();


    return levelSet;

}


int main(int argc, char* argv[]) {

    int numThreads = 1;

    if(argc != 1){
        numThreads = std::stoi(argv[1]);
    }

    std::cout << "running program with " << numThreads << " threads" << std::endl;


    omp_set_num_threads(numThreads);

    int numRuns = 10;

    NumericType gridDelta = 0.05;

    NumericType radius = 100.;

    std::stringstream csvOutput;

    csvOutput << numRuns << std::endl;

    csvOutput << numThreads << std::endl;

    auto start = std::chrono::high_resolution_clock::now(); 

    auto stop = std::chrono::high_resolution_clock::now(); 

    std::vector<lsSmartPointer<lsDomain<double, D>>> levelSets;

    //___________________________START GEOMETRY INPUT

    //lsSmartPointer<lsDomain<double, D>> levelSet = makeSphere(gridDelta, radius);

    std::vector<NumericType> planeNormal = {0. , 0. , 1.};

    lsSmartPointer<lsDomain<double, D>> levelSet = makeTrench(gridDelta, planeNormal);

    //lsSmartPointer<lsDomain<double, D>> levelSet;
    
    //lsReader<NumericType, D>(levelSet, "../data/rawLS4-14.lvst").apply();

    levelSets.push_back(levelSet);  

    //lsMesh mesh;
    //std::cout << "Extracting surface mesh..." << std::endl;
    //lsToSurfaceMesh<NumericType, D>(*(levelSets.back()), mesh).apply();
    //lsToDiskMesh<NumericType, D>(*(levelSets.back()), mesh).apply();
    //lsVTKWriter(mesh, lsFileFormatEnum::VTU , "mesh").apply();

    //___________________________END GEOMETRY INPUT

    std::cout << "num segments: " << levelSets.back()->getNumberOfSegments() << std::endl;

    std::cout << "FMM..." << std::endl;

    //lsEikonalExpandTest<NumericType, D> expander(levelSets.back());

    lsExpand<NumericType, D> expander(levelSets.back(), 3);

    expander.apply(); 

    double timeSum = 0.;

    std::cout << "num segments: " << levelSets.back()->getNumberOfSegments() << std::endl;

    std::cout << "Flagging Tests" << std::endl; 

    myFlagging<NumericType, D> myFlagger(levelSets.back(), 1e-3);

    std::cout << "Flagging Tests Shape Operator" << std::endl;

    csvOutput << "Shape Operator" << std::endl;

    for(int i =0; i < numRuns; i++){

        //start = std::chrono::high_resolution_clock::now(); 

        double time = myFlagger.apply(0);

        //stop = std::chrono::high_resolution_clock::now(); 

        std::cout << time << " "; 

        csvOutput << time << ";";

        timeSum += time;

    }
    csvOutput << std::endl;
    std::cout << std::endl;

    std::cout << "My Flagging Shape AVG: " << timeSum/numRuns << std::endl; 

    myFlagger.createFlagsOutput(0);

    myFlagging<NumericType, D> myFlagger2(levelSets.back(), 1e-4);

    timeSum = 0;

    std::cout << "Flagging Tests General Formula" << std::endl;

    csvOutput << "General Formula" << std::endl;

    for(int i =0; i < numRuns; i++){

        //start = std::chrono::high_resolution_clock::now(); 

        double time = myFlagger2.apply(1);

        //stop = std::chrono::high_resolution_clock::now(); 

        std::cout << time << " "; 

        csvOutput << time << ";";

        timeSum += time;

    }
    csvOutput << std::endl;
    std::cout << std::endl;

    std::cout << "My Flagging General AVG: " << timeSum/numRuns << std::endl; 

    //myFlagger2.createFlagsOutput(1);

/*    std::cout << "Flagging Test Normals" << std::endl;

    csvOutput << "Normal Flagging" << std::endl;

    silvacoLikeFlagging<NumericType,D> silvacoFlagger(levelSets.back(), 0.16);

    for(int i =0; i < numRuns; i++){

        start = std::chrono::high_resolution_clock::now(); 

        silvacoFlagger.apply();

        stop = std::chrono::high_resolution_clock::now(); 

        std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << " "; 

        csvOutput << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << ";";

    }
    csvOutput << std::endl;
    std::cout << std::endl;
    
    silvacoFlagger.createFlagsOutput();

    std::cout << "Silvaco Flagging: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << std::endl; 
*/
    std::ofstream output;

    output.open("timings/timingsFlagging" + std::to_string(numThreads) + ".csv");

    output << csvOutput.rdbuf();

    output.close();


/*    auto narrowband3 = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(levelSet, narrowband3, true, true, 0.5).apply();

    //TODO: create output function!

    //narrowband.insertNextVectorData(normal, "Normal");
    //narrowband3.insertNextScalarData(curve, "curvature");
  
    lsVTKWriter(narrowband3, lsFileFormatEnum::VTU , "Flags" ).apply();

*/
    std::cout << "Finished" << std::endl;

    return 0;
}




