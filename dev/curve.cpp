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



//____________testing not necessary_________________

#include <chrono>
#include <lsCalculateNormalVectors.hpp>

//____________testing end___________________________




int main() {

  constexpr int D = 3;
  typedef double NumericType;

  double gridDelta = 0.5;

  omp_set_num_threads(1);

  /*double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[3];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  
  lsDomain<double,D> levelSet(bounds, boundaryCons, gridDelta);*/

  std::vector<lsDomain<double, D> *> levelSets;

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


  //create sphere
  
  lsDomain<double,D> levelSet(gridDelta);

  std::cout << "creating sphere..." << std::endl;
  double origin[3] = {5., 0., 0.};
  double radius = 7.0;

  lsMakeGeometry<double, D>(levelSet, lsSphere<double, D>(origin, radius)).apply();

  levelSets.push_back(&levelSet);  
  


  std::cout << "Number of points: " << levelSets.back()->getDomain().getNumberOfPoints() << std::endl;





  std::chrono::time_point<std::chrono::system_clock> start, end;

  //test



  //expand level set so that all points exist
  //test_curvature.prepareLS(*(levelSets.back()));

  
  //TODO: this is shit the ls has to be extended here so that the sparsebox iterator is defined
  lsCalculateCurvatures<double, D>::prepareLS(*(levelSets.back()));

  lsCalculateCurvatures<double, D> test_curvature(*(levelSets.back()));



  lsCalculateNormalVectors<double, D> test_normals(*(levelSets.back()));

  start = std::chrono::system_clock::now();
 
  test_normals.apply();

  end = std::chrono::system_clock::now();

  std::cout << "Time for normals: "<< std::endl << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;



  start = std::chrono::system_clock::now();
 
  test_curvature.apply();

  end = std::chrono::system_clock::now();

  std::cout << "Time for my iterator: "<< std::endl << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << std::endl;


  std::cout << "num of normals: " << levelSets.back()->getNormalVectors().size() << std::endl;


  // prepare curvature data to write into output file



  auto& gauss_curve = test_curvature.getGaussianCurvature();
  auto& mean_curve = test_curvature.getMeanCurvature();
  auto& my_normals = test_curvature.getNormals();
  auto& normalVectors = levelSets.back()->getNormalVectors();

  std::cout << "anz my " << my_normals.size() << std::endl;
  std::cout << "anz ot " << normalVectors.size() << std::endl;

  std::vector<double> gauss;
  std::vector<double> mean;
  std::vector<std::array<double, D>> my_normal;
  std::vector<std::array<double, D>> normal;





  for (hrleConstSparseIterator<lsDomain<double, D>::DomainType> it(
           levelSets.back()->getDomain());
       !it.isFinished(); ++it) {
    if (!it.isDefined() || std::abs(it.getValue()) > 0.5)
      continue;
    normal.push_back(normalVectors[it.getPointId()]);
    my_normal.push_back(my_normals[it.getPointId()]);
    gauss.push_back(gauss_curve[it.getPointId()]);
    mean.push_back(mean_curve[it.getPointId()]);


  }

  std::cout << "anz gauss after loop " << gauss.size() << std::endl;

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



  lsMesh pointcloud;
  std::cout << "Extracting point cloud..." << std::endl;
  lsToDiskMesh<double, D>(*(levelSets.back()), pointcloud).apply();

  pointcloud.insertNextScalarData(gauss , "gaussian curvature");
  pointcloud.insertNextScalarData(mean , "mean curvature");

  lsVTKWriter(pointcloud, lsFileFormatEnum::VTU , "sphere_points").apply();

  //std::cout << "num of curves: " << test_curvature.getGaussianCurvature().size() << std::endl;




  std::cout << "Finished" << std::endl;
 

  return 0;
}


 


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

