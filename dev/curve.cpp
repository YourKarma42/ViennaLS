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





int main() {

  constexpr int D = 3;
  omp_set_num_threads(1);


  double gridDelta = 0.5;


  //create sphere
  std::cout << "creating sphere..." << std::endl;
  double origin[3] = {5., 0., 0.};
  double radius = 7.3;
  lsDomain<double, D> sphere1(gridDelta); //, boundaryCons);
  lsMakeGeometry<double, D>(sphere1, lsSphere<double, D>(origin, radius))
      .apply();


  std::cout << "Number of points: " << sphere1.getDomain().getNumberOfPoints()
           << std::endl;

  
  int num_points = sphere1.getDomain().getNumberOfPoints();

  std::vector<double> test(num_points, 13.37);


  lsMesh mesh;
  std::cout << "Extracting..." << std::endl;
  lsToSurfaceMesh<double, D>(sphere1, mesh).apply();

  mesh.insertNextScalarData(test, "test_data");

  lsVTKWriter(mesh, lsFileFormatEnum::VTU , "sphere").apply();

  lsMesh narrowband;
  std::cout << "Extracting..." << std::endl;
  lsToMesh<double, D>(sphere1, narrowband).apply();

  narrowband.insertNextScalarData(test, "test_data");

  lsVTKWriter(narrowband, lsFileFormatEnum::VTU ,"sphere_narrowband").apply();

  lsMesh pointcloud;
  std::cout << "Extracting..." << std::endl;
  lsToDiskMesh<double, D>(sphere1, pointcloud, 10.0).apply();

  pointcloud.insertNextScalarData(test, "test_data");

  lsVTKWriter(pointcloud, lsFileFormatEnum::VTU ,"sphere_points").apply();


  std::cout << "Finished" << std::endl;

 

  return 0;
}
