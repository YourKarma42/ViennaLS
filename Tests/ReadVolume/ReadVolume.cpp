#include <algorithm>
#include <iostream>

#include <lsDomain.hpp>
#include <lsFromVolumeMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

int main() {
  constexpr int D = 2;
   

  double gridDelta = 5e-10;
  double bounds[2 * D] = {-3.5e-8, 3.5e-8, -5e-8, 5e-8};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  boundaryCons[1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  std::vector<lsDomain<double, D>> levelSets(
      5, lsDomain<double, D>(bounds, boundaryCons, gridDelta));

  // Read mesh
  lsMesh initialMesh;
  lsVTKReader(initialMesh, "initial.vtk").apply();
  initialMesh.print();

  lsFromVolumeMesh<double, D>(levelSets, initialMesh).apply();

  for (unsigned i = 0; i < levelSets.size(); ++i) {
    lsMesh mesh;
    lsToSurfaceMesh<double, D>(levelSets[i], mesh).apply();
    lsVTKWriter(mesh, "LSsurface-" + std::to_string(i) + ".vtk").apply();
  }

  return 0;
}
