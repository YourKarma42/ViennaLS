#ifndef LS_TO_DISK_MESH_HPP
#define LS_TO_DISK_MESH_HPP

#include <lsPreCompileMacros.hpp>

#include <hrleSparseIterator.hpp>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMesh.hpp>

/// This class creates a mesh from the level set
/// with all grid points with a level set value <= 0.5.
/// These grid points are shifted in space towards the
/// direction of their normal vector by grid delta * LS value.
/// Grid delta and the origin grid point are saved for each point.
/// This allows for a simple setup of disks for ray tracing.
template <class T, int D> class lsToDiskMesh {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  lsDomain<T, D> *levelSet = nullptr;
  lsMesh *mesh = nullptr;
  T maxValue = 0.5;

public:
  lsToDiskMesh() {}

  lsToDiskMesh(lsDomain<T, D> &passedLevelSet, lsMesh &passedMesh,
               T passedMaxValue = 0.5)
      : levelSet(&passedLevelSet), mesh(&passedMesh), maxValue(passedMaxValue) {
  }

  void setLevelSet(lsDomain<T, D> &passedLevelSet) {
    levelSet = &passedLevelSet;
  }

  void setMesh(lsMesh &passedMesh) { mesh = &passedMesh; }

  void setMaxValue(const T passedMaxValue) { maxValue = passedMaxValue; }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsToDiskMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsToDiskMesh.")
          .print();
      return;
    }

    mesh->clear();

    maxValue = std::abs(maxValue);

    lsExpand<T, D>(*levelSet, (maxValue * 4) + 1).apply();
    lsCalculateNormalVectors<T, D>(*levelSet, maxValue, false).apply();
    const auto &normalVectors = levelSet->getNormalVectors();

    const T gridDelta = levelSet->getGrid().getGridDelta();

    // set up data arrays
    std::vector<double> values;
    std::vector<double> gridSpacing;
    std::vector<std::array<double, 3>> normals;

    for (hrleConstSparseIterator<hrleDomainType> it(levelSet->getDomain());
         !it.isFinished(); ++it) {

      if (!it.isDefined() || std::abs(it.getValue()) > maxValue) {
        continue;
      }

      // insert corresponding node shifted by ls value in direction of the
      // normal vector
      std::array<double, 3> normal = {0, 0, 0};
      for(unsigned i = 0; i < D; ++i) {
        normal[i] = normalVectors[it.getPointId()][i];
      }
      double modulus = 0.;
      for (unsigned i = 0; i < D; ++i) {
        modulus += normal[i] * normal[i];
      }
      modulus = std::sqrt(modulus);
      for (auto &n : normal) {
        n /= modulus;
      }
      normals.push_back(normal);

      std::array<double, 3> node;
      node[2] = 0.;
      for (unsigned i = 0; i < D; ++i) {
        // original position
        node[i] = double(it.getStartIndices(i)) * gridDelta;
        // shift to surface
        node[i] -= it.getValue() * gridDelta * normal[i] / modulus;
      }

      // insert vertex
      std::array<unsigned, 1> vertex;
      vertex[0] = mesh->nodes.size();
      mesh->insertNextVertex(vertex);

      mesh->insertNextNode(node);

      // add data into mesh
      values.push_back(it.getValue());
      gridSpacing.push_back(gridDelta);
    }

    mesh->insertNextScalarData(values, "LSValues");
    mesh->insertNextScalarData(gridSpacing, "gridSpacing");
    mesh->insertNextVectorData(normals, "Normals");
  }
};

#endif // LS_TO_DISK_MESH_HPP
