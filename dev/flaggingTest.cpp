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



#include <lsConvertEuclid.hpp>
#include <lsEikonalExpand.hpp>

#include <../dev/silvacoLikeFlagging.hpp>
#include <../dev/myFlagging.hpp>


#include <omp.h>



//____________testing not necessary_________________

#include <chrono>

#include <lsCalculateNormalVectors.hpp>

//____________testing end___________________________

constexpr int D = 3;
typedef double NumericType;

lsDomain<double, D> makeTrench(double gridDelta){

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

    lsDomain<double, D> levelSet(bounds, boundaryCons, gridDelta);

    double origin[3] = {0., 0., 0.};
    double planeNormal[3] = {0., D == 2, D == 3};

    //create the plane
    lsMakeGeometry<double, D>(levelSet, lsPlane<double, D>(origin, planeNormal))
        .apply();

    {
        // create layer used for booling
        std::cout << "Creating box..." << std::endl;
        lsDomain<double, D> trench(bounds, boundaryCons, gridDelta);
        double minCorner[3] = {-extent - 1, -extent / 4., -15.};
        double maxCorner[3] = {extent + 1, extent / 4., 1.0};
        lsMakeGeometry<double, D>(trench, lsBox<double, D>(minCorner, maxCorner))
            .apply();


        // Create trench geometry
        std::cout << "Booling trench..." << std::endl;
        lsBooleanOperation<double, D>(levelSet, trench,
                                    lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
            .apply();
    }

    return levelSet;
    
}

lsDomain<double, D> makeSphere(double gridDelta, double radius){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0., 0., 0.};
    
    lsDomain<double,D> levelSet(gridDelta);

    lsMakeGeometry<double, D>(levelSet, lsSphere<double, D>(origin, radius)).apply();


    return levelSet;

}


int main(int argc, char* argv[]) {

    int numThreads = 1;

    if(argc != 1){
        numThreads = std::stoi(argv[1]);
    }

    std::cout << "running program with " << numThreads << " threads" << std::endl;


    omp_set_num_threads(numThreads);

    NumericType gridDelta = 0.125;

    auto start = std::chrono::high_resolution_clock::now(); 

    auto stop = std::chrono::high_resolution_clock::now(); 

    std::vector<lsDomain<NumericType, D> *> levelSets;

    //lsDomain<NumericType,D> levelSet = makeSphere(gridDelta, 50.);

    lsDomain<NumericType,D> levelSet = makeTrench(gridDelta);

    levelSets.push_back(&levelSet);  

    //lsMesh mesh;
    //std::cout << "Extracting surface mesh..." << std::endl;
    //lsToSurfaceMesh<NumericType, D>(*(levelSets.back()), mesh).apply();
    //lsToDiskMesh<NumericType, D>(*(levelSets.back()), mesh).apply();
    //lsVTKWriter(mesh, lsFileFormatEnum::VTU , "mesh").apply();



    //convert level set
    std::cout << "Converting..." << std::endl;

    lsConvertEuclid<NumericType, D>  converter(*(levelSets.back()));

    converter.apply();

    //get the active grid points of the level set
    auto activePoints = converter.getActivePoints();

    std::cout << "FMM..." << std::endl;

    lsEikonalExpand<NumericType, D> expander(*(levelSets.back()), activePoints);

    expander.apply(); 

    std::cout << "Flagging Tests" << std::endl; 

    myFlagging<NumericType, D> myFlagger(*(levelSets.back()), 1e-3);


    for(int i =0; i < 5; i++){

        start = std::chrono::high_resolution_clock::now(); 

        myFlagger.apply(activePoints);

        stop = std::chrono::high_resolution_clock::now(); 

        std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << " "; 

    }
    std::cout << std::endl;

    myFlagger.createFlagsOutput(activePoints);



    std::cout << "My Flagging: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << std::endl; 

    silvacoLikeFlagging<NumericType,D> silvacoFlagger(*(levelSets.back()), 0.34906585039887);

    for(int i =0; i < 5; i++){

        start = std::chrono::high_resolution_clock::now(); 

        silvacoFlagger.apply(activePoints);

        stop = std::chrono::high_resolution_clock::now(); 

        std::cout << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << " "; 

    }

    silvacoFlagger.createFlagsOutput(activePoints);

    std::cout << "Silvaco Flagging: " << std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count() << std::endl; 




    lsMesh narrowband3;
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(levelSet, narrowband3, true, true).apply(activePoints);

    //TODO: create output function!

    //narrowband.insertNextVectorData(normal, "Normal");
    //narrowband3.insertNextScalarData(curve, "curvature");
  
    lsVTKWriter(narrowband3, lsFileFormatEnum::VTU , "narrowband" ).apply();


    std::cout << "Finished" << std::endl;

    return 0;
}




