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


//#include <hrleDomain.hpp>




//#include <hrleSparseBoxIterator.hpp>
//#include <hrleVectorType.hpp>

#include<lsCalculateCurvatures.hpp>

#include <lsConvertEuclid.hpp>
#include <lsEikonalExpand.hpp>



//____________testing not necessary_________________


#include <omp.h>
#include <chrono>
#include <lsCalculateNormalVectors.hpp>

//____________testing end___________________________

constexpr int D = 3;
typedef double NumericType;


void createPointCloudOutput(lsDomain<NumericType, D> &passedlsDomain, std::string outputName){

  std::cout << "Calculating Curvatures..." << std::endl;


  std::chrono::time_point<std::chrono::system_clock> start, end;

  lsConvertEuclid<NumericType, D>  converter(passedlsDomain);

  converter.apply();

  auto activePoints = converter.getActivePoints();

  lsEikonalExpand<NumericType, D> expander(passedlsDomain, activePoints);

  expander.apply(); 
/*
  lsExpand<NumericType, D> expander(passedlsDomain, 8);

  expander.apply();
*/

  

  lsCalculateCurvatures<double, D> test_curvature(passedlsDomain);

  //lsCalculateNormalVectors<double, D> test_normals(passedlsDomain);

  start = std::chrono::system_clock::now();
 
  //test_normals.apply();

  end = std::chrono::system_clock::now();

  std::cout << "Time for normals: "<< std::endl << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;



  start = std::chrono::system_clock::now();
 
  test_curvature.apply();

  end = std::chrono::system_clock::now();

  std::cout << "Time for my calculation: "<< std::endl << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;




  auto& gauss_curve = test_curvature.getGaussianCurvature();
  auto& mean_curve = test_curvature.getMeanCurvature();
  auto& normGrads = test_curvature.getNormGrad();
  auto& my_normals = test_curvature.getNormals();
  //auto& normalVectors = passedlsDomain.getNormalVectors();



  std::vector<double> gauss;
  std::vector<double> mean;
  std::vector<double> normGrad;
  std::vector<std::array<double, D>> my_normal;
  std::vector<std::array<double, D>> normal;


  // prepare curvature data to write into output file
  for (hrleConstSparseIterator<lsDomain<double, D>::DomainType> it(
           passedlsDomain.getDomain());
       !it.isFinished(); ++it) {
    if (!it.isDefined() || std::abs(it.getValue()) > 0.5)
      continue;

    //TODO: think of more elegant solution
    //if((activePoints.find(it.getStartIndices()) != activePoints.end())){

      gauss.push_back(std::abs(gauss_curve[it.getPointId()]));
      mean.push_back(std::abs(mean_curve[it.getPointId()]));
      normGrad.push_back(normGrads[it.getPointId()]);
      my_normal.push_back(my_normals[it.getPointId()]);      

    /*}else{

      gauss.push_back(0.);
      mean.push_back(0.);
      normGrad.push_back(0.);
      my_normal.push_back(my_normals[it.getPointId()]);
    }*/
    



  }





  std::cout << "Grid Delta: " << passedlsDomain.getGrid().getGridDelta() << std::endl;

  double sum = 0;

  for(auto c: gauss)
    sum += c;

  std::cout << "Gauss average " << (sum/gauss.size()) << std::endl;

  sum = 0;

  for(auto c: mean)
    sum += c;

  std::cout << "mean average " << (sum/gauss.size()) << std::endl;


  lsMesh pointcloud;
  std::cout << "Extracting point cloud..." << std::endl;
  lsToDiskMesh<double, D>(passedlsDomain, pointcloud).apply();

  pointcloud.insertNextScalarData(gauss , "Gaussian Curvature");
  pointcloud.insertNextScalarData(mean , "Mean Curvature");
  pointcloud.insertNextScalarData(normGrad , "Norm of Gradient");

  pointcloud.insertNextVectorData(my_normal, "my Normals");

  lsVTKWriter(pointcloud, lsFileFormatEnum::VTU , "points" + outputName ).apply();



}

NumericType testF(NumericType x, NumericType y){
  std::cout << "testitest in function" << std::endl;
  return 0;
}


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

lsDomain<double, D> makeSphere(double gridDelta){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0., 0., 0.};
    double radius = 10.0;

    lsDomain<double,D> levelSet(gridDelta);

    lsMakeGeometry<double, D>(levelSet, lsSphere<double, D>(origin, radius)).apply();

    return levelSet;

}



int main() {


  omp_set_num_threads(1);

  //double gridDelta = 0.5;

  std::vector<double> gridDeltas = { 0.5, 0.25, 0.125};

  std::vector<std::vector<NumericType>> additionalParams;
  /*
  //90 grad
  std::vector<NumericType> planeNormal = {0. , 1. , 0.};
  additionalParams.push_back(planeNormal);

  //22.5 grad
  planeNormal = {0.0, 0.38268343, 0.92387953};
  additionalParams.push_back(planeNormal);

  //45 grad
  planeNormal = {0.0, 0.70710678, 0.70710678};
  additionalParams.push_back(planeNormal);

  //67.5 grad
  planeNormal ={0.0, 0.92387953, 0.38268343};
  additionalParams.push_back(planeNormal);
  */


  /*
  std::vector<lsDomain<double, D> *> levelSets;
    
  lsDomain<double,D> levelSet(0.5);

  lsMakeGeometry<double, D>(levelSet, 
    lsGeometryByFunction<NumericType, D>(std::function <NumericType (NumericType, NumericType)>(testF)) ).apply();

  int a;

  std::cin >> a;
  */


  for(auto gridDelta: gridDeltas){

    int i=0;

    //for(auto param: additionalParams){

  double gridDelta = 0.5;


      //lsDomain<double,D> levelSet = makePlane(gridDelta, param);

     // lsDomain<double,D> levelSet = makeSphere(gridDelta);

      //levelSets.push_back(&levelSet);  
      

      //createPointCloudOutput(*(levelSets.back()), "sphere" + std::to_string(gridDelta) + " " + std::to_string(i));

      i++;
    //}

  }

  std::cout << "Finished" << std::endl;
 

  return 0;
}

  //create trench/plane

  /*
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
  }*/

 
/*
  lsMesh narrowband;
  std::cout << "Extracting narrowband..." << std::endl;
  lsToMesh<double, D>(*(levelSets.back()), narrowband, true, true).apply();

  narrowband.insertNextScalarData(gauss , "gaussian curvature");

  narrowband.insertNextScalarData(mean , "mean curvature");

  narrowband.insertNextVectorData(my_normal, "my Normals");

  narrowband.insertNextVectorData(normal, "Normals");

  lsVTKWriter(narrowband, lsFileFormatEnum::VTU ,"sphere_narrowband").apply();

  lsMesh mesh;
  std::cout << "Extracting surface mesh..." << std::endl;
  lsToSurfaceMesh<double, D>(*(levelSets.back()), mesh).apply();
  //mesh.insertNextScalarData(gauss , "gaussian curvature");
  lsVTKWriter(mesh, lsFileFormatEnum::VTU , "sphere_mesh").apply();
*/


  //output save for later
  /*int num_points = sphere1.getDomain().getNumberOfPoints();

  std::vector<double> test(num_points, 13.37);


  lsMesh mesh;
  std::cout << "Extracting surface mesh..." << std::endl;
  lsToSurfaceMesh<double, D>(sphere1, mesh).apply();

  //mesh.insertNextScalarData(test, "test_data");

  lsVTKWriter(mesh, "sphere").apply();

  lsVTKWriter(mesh, lsFileFormatEnum::VTU ,"sphere").apply();

  lsMesh narrowband;
  std::cout << "Extracting narrowband..." << std::endl;
  lsToMesh<double, D>(sphere1, narrowband).apply();

  //narrowband.insertNextScalarData(test, "test_data");

  lsVTKWriter(narrowband, lsFileFormatEnum::VTU ,"sphere_narrowband").apply();

  lsMesh pointcloud;
  std::cout << "Extracting point cloud..." << std::endl;
  lsToDiskMesh<double, D>(sphere1, pointcloud, 1.0).apply();

  pointcloud.insertNextScalarData(test, "test_data");

  lsVTKWriter(pointcloud, lsFileFormatEnum::VTU ,"sphere_points").apply();*/

