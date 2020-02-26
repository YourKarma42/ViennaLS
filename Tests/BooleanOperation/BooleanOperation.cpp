#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Example of boolean operations on level sets
  using two spheres.
  \example BooleanOperation.cpp
*/

int main() {

  constexpr int D = 3;

  double gridDelta = 1.0;
  double bounds[2 * D] = {-20, 20, -20, 20};
  if (D == 3) {
    bounds[4] = -20;
    bounds[5] = 20;
  }
  typename lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  lsDomain<double, D> sphere1(bounds, boundaryCons, gridDelta);
  lsDomain<double, D> sphere2(bounds, boundaryCons, gridDelta);

  double origin[3] = {0., 0., 0.};
  double radius = 15.7;

  lsMakeGeometry<double, D>(sphere1, lsSphere<double, D>(origin, radius))
      .apply();
  origin[0] = 15.0;
  radius = 9.5;
  lsMakeGeometry<double, D>(sphere2, lsSphere<double, D>(origin, radius))
      .apply();

  {
    lsMesh mesh1, mesh2;

    std::cout << "Extracting..." << std::endl;
    lsToSurfaceMesh<double, D>(sphere1, mesh1).apply();
    lsToSurfaceMesh<double, D>(sphere2, mesh2).apply();

    lsVTKWriter(mesh1, "sphere1.vtk").apply();
    lsVTKWriter(mesh2, "sphere2.vtk").apply();

    lsToMesh<double, D>(sphere1, mesh1).apply();
    lsToMesh<double, D>(sphere2, mesh2).apply();

    lsVTKWriter(mesh1, "LS1.vtk").apply();
    lsVTKWriter(mesh2, "LS2.vtk").apply();
  }

  // Perform a boolean operation
  lsBooleanOperation<double, D>(sphere1, sphere2, lsBooleanOperationEnum::UNION)
      .apply();

  std::cout << "Extracting..." << std::endl;
  lsMesh mesh;
  lsToSurfaceMesh<double, D>(sphere1, mesh).apply();

  mesh.print();

  lsVTKWriter(mesh, "after.vtk").apply();

  return 0;
}
