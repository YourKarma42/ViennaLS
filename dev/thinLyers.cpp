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


#include <omp.h>


//____________testing not necessary_________________

#include <chrono>


//#include <lsAdvect.hpp>

//____________testing end___________________________

constexpr int D = 2;
typedef double NumericType;
typedef typename lsDomain<NumericType, D>::DomainType hrleDomainType;


void printLayers(std::vector<lsSmartPointer<lsDomain<double, D>>> layers){

    std::cout << "Extracting Layers as Meshes..." << std::endl;

    int i = 1;

    for(auto layer: layers){

        auto mesh = lsSmartPointer<lsMesh>::New();
        lsToSurfaceMesh<double, D>(layer, mesh).apply();
        lsVTKWriter(mesh, "Layer" + std::to_string(i)  + ".vtk").apply();

        i++;

    }



}


lsSmartPointer<lsDomain<double, D>> makeMaterialLayer(double gridDelta, double origin[3], double extent,
                            lsNormalizations normalization){

    std::cout << "Creating Material Layer..." << std::endl;

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
        lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta, normalization);

    double planeNormal[3] = {0., 1., 0.};


    auto lsWithPlane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
    lsMakeGeometry<double, D>(levelSet, lsWithPlane).apply();



    return levelSet;

}

lsSmartPointer<lsDomain<double, D>> makeMaterialLayerBox(double gridDelta, double minCorner[3], double maxCorner[3], double extent,
                            lsNormalizations normalization){

    std::cout << "Creating Material Layer..." << std::endl;

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

    auto levelSet = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons,
                                                           gridDelta);

    auto box = lsSmartPointer<lsBox<double, D>>::New(minCorner, maxCorner);
    lsMakeGeometry<double, D>(levelSet, box).apply();



    return levelSet;

}

lsSmartPointer<lsDomain<double, D>> makeBox(double gridDelta, double minCorner[3], double extent){

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

    auto levelSet = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons,
                                                           gridDelta);

    double maxCorner[D] = {extent + 1 , extent + 1};
    
    if(D == 3){
        
        maxCorner[2] = extent;
    }

    auto box = lsSmartPointer<lsBox<double, D>>::New(minCorner, maxCorner);
    lsMakeGeometry<double, D>(levelSet, box).apply();

    return levelSet;
}


//Idea of bool operation

//This algorithm does not require redistancing it only iterates over ls values!

//expandLS because of manhatten

//LS A / LS B

//mybool (only relative complement)
    //keeps Ls values up to a certain value

// when the value of LS B is choosen this cell is "marked"

// "expand" "marked" cells in the direction of LS A

// mark LS A as effected by boolean operation



int main() {


    omp_set_num_threads(1);

    double gridDelta = 0.0625;

    double extent = 30;

    std::cout << "Creating Material Layers..." << std::endl;


    std::vector<lsSmartPointer<lsDomain<double, D>>> layers;

    std::vector<std::pair<NumericType, NumericType>> matPos;

    //first min; second max
    matPos.push_back(std::make_pair(-1., 0.)); 

    matPos.push_back(std::make_pair(-3., -1.));

    matPos.push_back(std::make_pair(-4.5, -3.)); 

    matPos.push_back(std::make_pair(-5.2, -4.5));  

    matPos.push_back(std::make_pair(-5.8, -5.2)); 

    matPos.push_back(std::make_pair(-10., -5.8)); 

    double min[3] = {-extent - 1, 0., 0.};

    double max[3] = {extent + 1, 0., 0.};

    double posSingleLayer[3] = {0., 0., 0.};

    int maxNumLayers = 1;

    for(auto pos: matPos){
    
        //min[1] = pos.first;

        //max[1] = pos.second;
    
        //layers.push_back(makeMaterialLayerBox(gridDelta, min, max, extent, lsNormalizations::MANHATTEN));

        posSingleLayer[1] = pos.second;

        layers.push_back(makeMaterialLayer(gridDelta, posSingleLayer , extent, lsNormalizations::MANHATTEN));
    
        maxNumLayers++;
    }


    std::cout << "Booling Layers..." << std::endl;

    double minCorner[3] = {0., 0., 0.};

    minCorner[1] = matPos[matPos.size()-2].first;

    auto box = makeBox(gridDelta, minCorner, extent);

    for(auto layer: layers){

        lsBooleanOperation<double, D>(layer, box, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
            .apply();

    }


    printLayers(layers);

    std::cout << "Finished" << std::endl;

    return 0;
}


/*
    OLD BOOL LAYER

    std::cout << "Booling Layers..." << std::endl;

    double minCorner[3] = {0., 0., 0.};

    int numLayers = 1;

    for(auto pos: matPos){

        if(numLayers == maxNumLayers-1)
            break;

        minCorner[1] = pos.first;

        auto box = makeBox(gridDelta, minCorner, extent);

        //bool all previous layers
        for(int i=0; i < numLayers; i++){
            lsBooleanOperation<double, D>(layers[i], box, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
                .apply();

        }

        numLayers++;

    }





*/




/*
    double origin[3] = {0.0, 0.0, 0.0};

    layers.push_back(makeMaterialLayer(gridDelta, origin, extent, lsNormalizations::MANHATTEN));

    origin[1] = -1.0;

    layers.push_back(makeMaterialLayer(gridDelta, origin, extent, lsNormalizations::MANHATTEN));

    origin[1] = -3.0;

    layers.push_back(makeMaterialLayer(gridDelta, origin, extent, lsNormalizations::MANHATTEN));
*/





