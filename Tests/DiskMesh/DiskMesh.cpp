#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to write a level set to an efficient disk mesh
  used for ray tracing.
  \example DiskMesh.cpp
*/

int main() {
  constexpr int D = 2;
  double gridDelta = 0.7;

  lsDomain<double, D> levelSet(gridDelta);

  const double radius = 27.3;
  const hrleVectorType<double, D> centre(5.2, 0.);

  lsMakeGeometry<double, 2>(levelSet, lsSphere<double, D>(centre, radius))
      .apply();

  lsDomain<double, D> box(gridDelta);
  double min[D] = {10.3, -5.1};
  double max[D] = {40.7, 5.4};
  lsMakeGeometry<double, 2>(box, lsBox<double, D>(min, max)).apply();

  lsBooleanOperation<double, 2>(levelSet, box, lsBooleanOperationEnum::UNION)
      .apply();

  lsMesh mesh;
  lsToDiskMesh<double, D>(levelSet, mesh).apply();

  mesh.print();

  lsVTKWriter(mesh, lsFileFormatEnum::VTP, "diskMesh.vtp").apply();
  lsVTKWriter(mesh, "diskMesh.vtk").apply();

  lsToSurfaceMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "surfaceMesh.vtk").apply();

  return 0;
}
