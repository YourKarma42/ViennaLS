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


        auto mesh1 = lsSmartPointer<lsMesh>::New();
        lsToMesh<double, D>(layer, mesh, true, false, 3.).apply();
        lsVTKWriter(mesh, "pointsLayer" + std::to_string(i)  + ".vtk").apply();

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

//only for development!
typedef NumericType T;

  void myTestBool(lsSmartPointer<lsDomain<T, D>> levelSetA, lsSmartPointer<lsDomain<T, D>> levelSetB) {
    auto &grid = levelSetA->getGrid();
    auto newlsDomain = lsSmartPointer<lsDomain<T, D>>::New(grid);
    typename lsDomain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
    typename lsDomain<T, D>::DomainType &domain = levelSetA->getDomain();

    newDomain.initialize(domain.getNewSegmentation(), domain.getAllocation());

    std::vector<std::vector<double>> markedCellsSegments(
    newlsDomain->getNumberOfSegments());
    double pointsPerSegment =
        double(2 * newlsDomain->getDomain().getNumberOfPoints()) /
        double(newlsDomain->getLevelSetWidth());

#pragma omp parallel num_threads(newDomain.getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      auto &domainSegment = newDomain.getDomainSegment(p);

      auto &markedCellsSegmentVector = markedCellsSegments[p];
      markedCellsSegmentVector.reserve(pointsPerSegment);

      hrleVectorType<hrleIndexType, D> currentVector =
          (p == 0) ? grid.getMinGridPoint()
                   : newDomain.getSegmentation()[p - 1];

      hrleVectorType<hrleIndexType, D> endVector =
          (p != static_cast<int>(newDomain.getNumberOfSegments() - 1))
              ? newDomain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      hrleConstSparseIterator<hrleDomainType> itA(levelSetA->getDomain(),
                                                  currentVector);
      hrleConstSparseIterator<hrleDomainType> itB(levelSetB->getDomain(),
                                                  currentVector);

      while (currentVector < endVector) {

        //(const T &currentValue = comp(itA.getValue(), itB.getValue());
//___________________________My Compare Test_______________________________________________
        T currentValue = 0.;
        bool mark = false;

        T A = itA.getValue();

        T B = itB.getValue();

        if(A >= -B){
            currentValue = A;
        }else{
            currentValue = -B;
            mark = true;
        }

//___________________________My Compare Test End____________________________________________
        if (currentValue != lsDomain<T, D>::NEG_VALUE &&
            currentValue != lsDomain<T, D>::POS_VALUE) {
          domainSegment.insertNextDefinedPoint(currentVector, currentValue);
          if(mark){
            markedCellsSegmentVector.push_back(1);
          }else{
            markedCellsSegmentVector.push_back(0);
          }
        } else {
          domainSegment.insertNextUndefinedPoint(
              currentVector, (currentValue < 0) ? lsDomain<T, D>::NEG_VALUE
                                                : lsDomain<T, D>::POS_VALUE);
        }

        switch (compare(itA.getEndIndices(), itB.getEndIndices())) {
        case -1:
          itA.next();
          break;
        case 0:
          itA.next();
        default:
          itB.next();
        }
        currentVector = std::max(itA.getStartIndices(), itB.getStartIndices());
      }
    }

    unsigned numberOfMarks = 0;
    for (unsigned i = 0; i < newlsDomain->getNumberOfSegments(); ++i) {
      numberOfMarks += markedCellsSegments[i].size();
    }
    markedCellsSegments[0].reserve(numberOfMarks);

    for (unsigned i = 1; i < newlsDomain->getNumberOfSegments(); ++i) {
      markedCellsSegments[0].insert(markedCellsSegments[0].end(),
                                    markedCellsSegments[i].begin(),
                                    markedCellsSegments[i].end());
    }
    // insert into pointData of levelSet
    auto &pointData = newlsDomain->getPointData();
    auto scalarDataPointer = pointData.getScalarData("BoolMarks");
    // if it does not exist, insert new marked cells
    if (scalarDataPointer == nullptr) {
      pointData.insertNextScalarData(markedCellsSegments[0], "BoolMarks");
    } else {
      // if it does exist, just swap the old with the new values
      *scalarDataPointer = std::move(markedCellsSegments[0]);
    }

    newDomain.finalize();
    newDomain.segment();
    //TODO: rethink this
    newlsDomain->setLevelSetWidth(levelSetA->getLevelSetWidth());
    
    //lsPrune<T, D>(newlsDomain).apply();

    levelSetA->deepCopy(newlsDomain);
  }






int main() {


    omp_set_num_threads(1);

    double gridDelta = 0.05;

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

        auto newLayer = makeMaterialLayer(gridDelta, posSingleLayer, extent, lsNormalizations::MANHATTEN);

        lsExpand(newLayer, 4.);

        layers.push_back(newLayer);
    
        maxNumLayers++;
    }


    std::cout << "Booling Layers..." << std::endl;

    double minCorner[3] = {0., 0., 0.};

    minCorner[1] = matPos[matPos.size()-2].first;

    auto box = makeBox(gridDelta, minCorner, extent);

    lsExpand(box, 4.);

    for(auto layer: layers){

        myTestBool(layer, box);

        //lsBooleanOperation<double, D>(layer, box, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        //    .apply();

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





