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


#include <lsWriter.hpp>
#include <lsReader.hpp>


#include <omp.h>



//____________testing not necessary_________________

#include <chrono>

#include <lsCalculateNormalVectors.hpp>

//____________testing end___________________________

constexpr int D = 3;
typedef double NumericType;

lsSmartPointer<lsDomain<double, D>> makeTrench(double gridDelta, std::vector<NumericType>& planeNormal){

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

    auto levelSet =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

    std::vector<NumericType> origin = {0., 0., 0.};

    {
      auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
      lsMakeGeometry<double, D>(levelSet, plane).apply();
    }

    {
        // create layer used for booling
        std::cout << "Creating box..." << std::endl;

        auto trench = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

        if(D == 3){
          double minCorner[3] = {-extent - 1, -extent / 4., -25.};
          double maxCorner[3] = {extent + 1, extent / 4., 1.0};
          auto box = lsSmartPointer<lsBox<double, D>>::New(minCorner, maxCorner);
          lsMakeGeometry<double, D>(trench, box).apply();
        }else{
          double minCorner[2] = {-extent / 4., -25.};
          double maxCorner[2] = {extent / 4., 1.0};
          auto box = lsSmartPointer<lsBox<double, D>>::New(minCorner, maxCorner);
          lsMakeGeometry<double, D>(trench, box).apply();
        }
        //lsMakeGeometry<double, D>(trench, lsBox<double, D>(minCorner, maxCorner))
        //    .apply();


        // Create trench geometry
        std::cout << "Booling trench..." << std::endl;
        lsBooleanOperation<double, D>(levelSet, trench,
                                      lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
    }

    return levelSet;
    
}

lsSmartPointer<lsDomain<double, D>> makeSphere(double gridDelta, double radius){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0.0, 0.0, 0.0};
    
    auto levelSet =
        lsSmartPointer<lsDomain<double, D>>::New(gridDelta);

    lsMakeGeometry<double, D>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();


    return levelSet;

}



int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();

  const double radius = 7.3;
  const hrleVectorType<double, D> centre(5., 0.);

  lsMakeGeometry<double, D>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(centre, radius))
      .apply();
/*
  lsPointData &data = levelSet->getPointData();
  typename lsPointData::ScalarDataType scalars;
  typename lsPointData::VectorDataType vectors;
  for (unsigned i = 0; i < levelSet->getNumberOfPoints(); ++i) {
    scalars.push_back(i);
    vectors.push_back(
        typename lsPointData::VectorDataType::value_type({double(i)}));
  }

  data.insertNextScalarData(scalars, "myScalars");
  data.insertNextVectorData(vectors, "llaalalalalaalalalalalaalal");
*/
  lsWriter<double, D>(levelSet, "test.lvst").apply();

  // read it in again
  auto newLevelSet = lsSmartPointer<lsDomain<double, D>>::New();
  lsReader<double, D>(newLevelSet, "test.lvst").apply();

  auto mesh = lsSmartPointer<lsMesh>::New();
  lsToSurfaceMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "test.vtk").apply();

  auto narrowband = lsSmartPointer<lsMesh>::New();
  std::cout << "Extracting narrowband..." << std::endl;
  lsToMesh<NumericType, D>(newLevelSet, narrowband, true, true, 0.5).apply();

  lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "loadedLvlSet" ).apply();

  return 0;
}



