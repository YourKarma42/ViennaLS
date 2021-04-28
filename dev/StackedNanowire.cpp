#include <iostream>
#include <array>
#include <deque>
#include <fstream>

#include <lsSmartPointer.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToDiskMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsAdvect.hpp>

#include <lsWriter.hpp>

#include <lsCalculateNormalVectors.hpp>

#include "derivatives.hpp"

constexpr int D = 3;
using NumericType = double;
constexpr NumericType gridDelta = 0.49999999;
unsigned outputNum = 0;

typedef typename lsDomain<NumericType, D>::DomainType hrleDomainType;

bool exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

// Implement own velocity field
class isotropic : public lsVelocityField<NumericType> {
  std::vector<double> velocities;
public:
  isotropic(std::vector<double> vel) : velocities(vel) {};

  double getScalarVelocity(
      const std::array<NumericType, 3> & /*coordinate*/, int material,
      const std::array<NumericType, 3> &/*normalVector = hrleVectorType<NumericType, 3>(0.)*/) {
    return (material < int(velocities.size()))?velocities[material]:0;
  }
};

// Directional etch for one material
class directional : public lsVelocityField<NumericType> {
  const std::array<NumericType, D> direction;
  std::vector<double> velocities;
  const double isoVelocity;
public:
  directional(std::array<NumericType, D> dir, std::vector<double> vel, double isoVel = 0) : direction(dir), velocities(vel), isoVelocity(isoVel) {};

  std::array<double, D> getVectorVelocity(
      const std::array<NumericType, 3> & coordinate, int material,
      const std::array<NumericType, 3> & normalVector = std::array<NumericType, 3>({})) {
    if(material < int(velocities.size())) {
      std::array<NumericType, D> dir(direction);
      for(unsigned i = 0; i < D; ++i) {
        dir[i] *= velocities[material];
        dir[i] += isoVelocity;
      }
      //if(coordinate == std::array<NumericType, 3>({10, 0, 28})){
      //  std::cout << "OUTPUT " << outputNum << " -- > Material: " << material << std::endl;
      //  std::cout << "Normal: " << normalVector[0] << ", " << normalVector[1] << ", " << normalVector[2] << std::endl;
      //  std::cout << "Velocity: " << dir[0] << ", " << dir[1] << ", " << dir[2] << std::endl;
      //}
      return dir;
    } else {
      return {0};
    }
  }
};

std::vector<std::vector<NumericType>> calcCurve(lsSmartPointer<lsDomain<NumericType, D>>& domain){

    hrleSparseBoxIterator<hrleDomain<NumericType, D>> neighborIterator(domain->getDomain(), 1);

    std::vector<std::vector<NumericType>> curvaturesReturn;
    std::vector<NumericType> curvatures;
    std::vector<NumericType> curvaturesBias;
    std::vector<NumericType> curvaturesAbs;

    //std::cout << "blub";

    lsExpand<NumericType, D>(domain, 7).apply();

    //std::cout << "blub2";

    NumericType gridDelta = domain->getGrid().getGridDelta();

    curvaturGeneralFormula<NumericType, D> generalFormula(gridDelta);

    curvaturGeneralFormulaBigStencilBias<NumericType, D> generalFormulaBias(gridDelta);

    for(hrleConstSparseIterator<hrleDomainType> centerIt(domain->getDomain());
      !centerIt.isFinished(); ++centerIt){

        if (!centerIt.isDefined() || (std::abs(centerIt.getValue()) > 0.5)) {
          continue;
        } 

        neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

        NumericType c = generalFormula(neighborIterator);

        curvatures.push_back(c); 

        curvaturesAbs.push_back(std::abs(c));

        curvaturesBias.push_back(generalFormulaBias(neighborIterator));   

    }

    curvaturesReturn.push_back(curvatures);
    curvaturesReturn.push_back(curvaturesBias);
    curvaturesReturn.push_back(curvaturesAbs);

    return curvaturesReturn;

}

void writeFeaturOutput(std::deque<lsSmartPointer<lsDomain<NumericType, D>>>& domains) {

  static unsigned numMat = 0;

  std::string path = "SN_Output/";

  for(unsigned i = 0; (i < domains.size()) || (i < numMat); ++i) {
    auto pointMesh = lsSmartPointer<lsMesh>::New();
    if(numMat != 0 && i >= numMat) {  // write emtpy meshes before so paraview displays it correctly
      for(unsigned j = 0; j < outputNum; ++j) {
        std::string fileName = std::to_string(i) + "-" + std::to_string(j) + ".vtk";
        if(!exists(fileName)){
          lsVTKWriter(pointMesh, path + "points-features" + fileName).apply();
        }
      }
    }
    
    if(i < domains.size()) {

      lsExpand<NumericType, D>(domains[i], 4).apply();

      //lsCalculateNormalVectors<NumericType, D>(domains[i]).apply();
      lsToDiskMesh<NumericType, D>(domains[i], pointMesh).apply();
      std::vector<std::vector<NumericType>> curves = calcCurve(domains[i]);       
      //lsToMesh<NumericType, D>(domains[i], pointMesh, true, true).apply();
      pointMesh->insertNextScalarData(curves[0], "curvature");
      pointMesh->insertNextScalarData(curves[1], "curvatureBias");
      pointMesh->insertNextScalarData(curves[2], "curvatureAbs");
    }

    std::string pointName = "points-features" + std::to_string(i) + "-" + std::to_string(outputNum) + ".vtk";
    lsVTKWriter(pointMesh, path + pointName).apply();
    
  }
  numMat = (domains.size()<numMat)?numMat:domains.size();
  // increase count
  //++outputNum;


}

void writeSurfaces(std::deque<lsSmartPointer<lsDomain<NumericType, D>>>& domains, bool writePoints = false) {

  static unsigned numMat = 0;

  std::string path = "SN_Output/";

  for(unsigned i = 0; (i < domains.size()) || (i < numMat); ++i) {
    auto mesh = lsSmartPointer<lsMesh>::New();
    auto pointMesh = lsSmartPointer<lsMesh>::New();
    if(numMat != 0 && i >= numMat) {  // write emtpy meshes before so paraview displays it correctly
      for(unsigned j = 0; j < outputNum; ++j) {
        std::string fileName = std::to_string(i) + "-" + std::to_string(j) + ".vtk";
        if(!exists(fileName)){
          lsVTKWriter(mesh, path + "surface-m" + fileName).apply();
          lsVTKWriter(mesh, path + "points-m" + fileName).apply();
        }
      }
    }
    
    if(i < domains.size()) {
      lsToSurfaceMesh<NumericType, D>(domains[i], mesh).apply();
      if(writePoints) {
        lsToMesh<NumericType, D>(domains[i], pointMesh).apply();
      }
    }
    std::string fileName = "surface-m" + std::to_string(i) + "-" + std::to_string(outputNum) + ".vtk";
    lsVTKWriter(mesh, path + fileName).apply();
    if(writePoints) {
      

      std::string pointName = "points-m" + std::to_string(i) + "-" + std::to_string(outputNum) + ".vtk";
      lsVTKWriter(pointMesh, path + pointName).apply();
    }

    //std::string levelSetName = "rawLS" + std::to_string(i) + "-" + std::to_string(outputNum) + ".lvst";
    //lsWriter<NumericType, D>(domains[i], path + levelSetName).apply();
  }

  writeFeaturOutput(domains);
  numMat = (domains.size()<numMat)?numMat:domains.size();
  // increase count
  ++outputNum;


}

class Process {
  public:
  std::string name;
  double time;
  lsSmartPointer<lsVelocityField<NumericType>> velocities;
  bool newLayer;

  template<class VelocityField>
  Process(std::string processName, double processTime, VelocityField processVelocities, bool newMaterial = false) : name(processName), time(processTime), newLayer(newMaterial) {
    velocities = std::dynamic_pointer_cast<lsVelocityField<NumericType>>(processVelocities);
  }
};

void execute(std::deque<lsSmartPointer<lsDomain<NumericType, D>>>& domains, std::vector<Process>& processes) {
  // Process loop
  for(auto& it : processes) {
    std::cout << "Running " << it.name << " for " << it.time << "s" << std::endl;
    if(it.newLayer) { // if new layer, copy last level set
      std::cout << "Adding new layer" << std::endl;
      domains.push_back(lsSmartPointer<lsDomain<NumericType, D>>::New(domains.back()));
    }

    auto advectionKernel = lsSmartPointer<lsAdvect<NumericType, D>>::New();
    advectionKernel->setVelocityField(it.velocities);
    for(auto &it : domains) {
      advectionKernel->insertNextLevelSet(it);
    }
      
    advectionKernel->setAdvectionTime(it.time);
    advectionKernel->apply();

    writeSurfaces(domains, true);
  }
}

void planarise(std::deque<lsSmartPointer<lsDomain<NumericType, D>>>& domains, double height) {
  auto plane = lsSmartPointer<lsDomain<NumericType, D>>::New(domains.front()); // copy domain
  double origin[3] = {0., (D==2)?height:0., (D==3)?height:0.};
  double planeNormal[3] = {0., (D==2)?-1:0., (D==3)?-1:0.};
  lsMakeGeometry<NumericType, D>(plane, lsSmartPointer<lsPlane<NumericType, D>>::New(origin, planeNormal)).apply();
  // now remove plane from domain
  for(auto &it : domains) {
    lsBooleanOperation<NumericType, D>(it, plane, lsBooleanOperationEnum::RELATIVE_COMPLEMENT).apply();
  }
  writeSurfaces(domains, true);
}

int main() {
  omp_set_num_threads(4);

  
  NumericType bounds[2 * D] = {0, 70, 0, 100, 0, 70}; // in nanometres
  lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = lsDomain<NumericType, D>::BoundaryType::PERIODIC_BOUNDARY;
  boundaryCons[1] = lsDomain<NumericType, D>::BoundaryType::PERIODIC_BOUNDARY;
  boundaryCons[2] = lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  std::deque<lsSmartPointer<lsDomain<NumericType, D>>> domains;
  domains.push_back(lsSmartPointer<lsDomain<NumericType, D>>::New(bounds, boundaryCons, gridDelta));

  { // Initialise domain with a single silicon plane ( at z= 10 because it is 10 nm high)
    NumericType origin[D] = {0., 0., 10.};
    NumericType planeNormal[D] = {0., 0., 1.};
    auto plane = lsSmartPointer<lsPlane<NumericType, D>>::New(origin, planeNormal);
    lsMakeGeometry<NumericType, D>(domains[0], plane).apply();
  }

  writeSurfaces(domains, true);

  // Set up first bunch of processes
  std::vector<Process> processes;

  auto isoVelocity = lsSmartPointer<isotropic>::New(std::vector<double>(10, 1.0));
  processes.push_back(Process("Si-Epitaxy", 7, isoVelocity, true));
  processes.push_back(Process("SiGe-Epitaxy", 8, isoVelocity, true));
  processes.push_back(Process("Si-Epitaxy", 7, isoVelocity, true));

  execute(domains, processes);

  // Add double patterning mask
  {
    domains.push_front(lsSmartPointer<lsDomain<NumericType, D>>::New(bounds, boundaryCons, gridDelta));
    hrleVectorType<double, D> min(25, -10, 31.9);
    hrleVectorType<double, D> max(45, 110, 60);
    lsMakeGeometry<NumericType, D> geometryFactory(domains.front(), lsSmartPointer<lsBox<NumericType, D>>::New(min, max));
    bool ignoreBNC = true;
    geometryFactory.setIgnoreBoundaryConditions(ignoreBNC); // so that the mask is not mirrored inside domain at bounds
    geometryFactory.apply();

    // wrap all other layers accordingly
    for(unsigned i = 4; i < domains.size(); ++i) {
      lsBooleanOperation<NumericType, D>(domains[i], domains.front(), lsBooleanOperationEnum::UNION).apply();
    }
  }
  writeSurfaces(domains, true);

  // Double patterning processes
  {
    processes.clear();
    processes.push_back(Process("DP-Depo", 15, isoVelocity, true));
    std::array<NumericType, D> direction = {0, 0, -1};
    auto directVelocity = lsSmartPointer<directional>::New(direction, std::vector<double>({0,0,0,0,0,1}));
    processes.push_back(Process("DP-Patterning", 20, directVelocity));
  }

  execute(domains, processes);
  
  // Remove mask with boolean op
  {
    for(unsigned i = 4; i < domains.size(); ++i) {
      lsBooleanOperation<NumericType, D>(domains[i], domains.front(), lsBooleanOperationEnum::RELATIVE_COMPLEMENT).apply();
    }
    domains.pop_front(); // now remove
  }
  writeSurfaces(domains, true);


  // pattern si/sige/si stack
  {
    processes.clear();
    std::array<NumericType, D> direction = {0, 0, -1};
    auto directVelocity = lsSmartPointer<directional>::New(direction, std::vector<double>({0.5,1,1,1}));
    processes.push_back(Process("Si/Ge/Si-Patterning", 30, directVelocity));
  }
  execute(domains, processes);

  // Remove DP mask
  domains.pop_back();
  writeSurfaces(domains, true);


//return 0;
  // deposit dummy gate material
  {
    processes.clear();
    processes.push_back(Process("DG-Deposit", 55, isoVelocity, true));
  }
  execute(domains, processes);

  // dummy gate CMP at 80nm height
  planarise(domains, 80.0);

  // dummy gate mask addition
  {
    domains.push_front(lsSmartPointer<lsDomain<NumericType, D>>::New(bounds, boundaryCons, gridDelta));
    hrleVectorType<double, D> min(-10, 30, 75);
    hrleVectorType<double, D> max(80, 70, 90);
    lsMakeGeometry<NumericType, D> geometryFactory(domains.front(), lsSmartPointer<lsBox<NumericType, D>>::New(min, max));
    bool ignoreBNC = true;
    geometryFactory.setIgnoreBoundaryConditions(ignoreBNC); // so that the mask is not mirrored inside domain at bounds
    geometryFactory.apply();

    // wrap all other layers accordingly
    for(unsigned i = 5; i < domains.size(); ++i) { // only bool with 5th layer
      lsBooleanOperation<NumericType, D>(domains[i], domains.front(), lsBooleanOperationEnum::UNION).apply();
    }
  }
  writeSurfaces(domains, true);

  // dummy gate patterning
  {
    processes.clear();
    std::array<NumericType, D> direction = {0, 0, -1};
    auto directVelocity = lsSmartPointer<directional>::New(direction, std::vector<double>({0, 0, 0, 0, 0, 1}));
    processes.push_back(Process("DG-Patterning", 90, directVelocity));
  }
  execute(domains, processes);

  // Remove mask
  {
    domains.pop_front(); // now remove
  }
  writeSurfaces(domains, true);

  // spacer deposition, spacer patterning, fin patterning, SD Epitaxy, dielectric deposition
  {
    processes.clear();
    processes.push_back(Process("Spacer-Deposition", 12, isoVelocity, true));
    std::array<NumericType, D> direction = {0, 0, -1};
    auto spacerPatterning = lsSmartPointer<directional>::New(direction, std::vector<double>({0.1, 0.1, 0.1, 0.1, 0.1, 1}));
    processes.push_back(Process("Spacer-Patterning", 40, spacerPatterning));

    auto finPatterning = lsSmartPointer<directional>::New(direction, std::vector<double>({-0.05, 1, 1, 1, 0.1, 0.1})); // 0.1 because there needs to be some etching on the top layer for fins to be etched properly (apparently some numeric issues)
    processes.push_back(Process("Fin-Patterning", 19, finPatterning));
    
    auto SDEpitaxy = lsSmartPointer<isotropic>::New(std::vector<double>({1.0, 1.0, 1.0, 1.0, 0., 0., 1.0}));
    processes.push_back(Process("SD-Epitaxy", 11, SDEpitaxy, true));

    processes.push_back(Process("Dielectric-Deposition", 35, isoVelocity, true));
  }
  execute(domains, processes);

  /* iretative fin patterning */
  // {
  //   std::array<NumericType, D> direction = {0, 0, -1};
  //   auto finPatterning = lsSmartPointer<directional>::New(direction, std::vector<double>({-0.05, 1, 1, 1, 0.1, 0.1})); // 0.1 because there needs to be some etching on the top layer for fins to be etched properly (apparently some numeric issues)
  //   // processes.push_back(Process("Fin-Patterning", 19, finPatterning));

  //   /* Iterative outputs */
  //   auto advectionKernel = lsSmartPointer<lsAdvect<NumericType, D>>::New();
  //   advectionKernel->setVelocityField(finPatterning);
  //   advectionKernel->setSaveAdvectionVelocities(true);
  //   for(auto &it : domains) {
  //     advectionKernel->insertNextLevelSet(it);
  //   }
  //   double advectionTime = 19.0;
  //   double currentTime = 0.0;
  //   while (currentTime < advectionTime) {
  //     advectionKernel->apply(); // no advection time set, so only advect once
  //     currentTime += advectionKernel->getAdvectedTime();
  //     writeSurfaces(domains, true);
  //   }
  // }

  // dummy gate and dielectric CMP at 72.5nm height
  planarise(domains, 72.5);

  // remove dummy gate
  {
    processes.clear();
    auto dummyGateEtch = lsSmartPointer<isotropic>::New(std::vector<double>({-0.05, 0., 0., 0., -1.}));
    processes.push_back(Process("DummyGate-Removal", 80, dummyGateEtch));
  }
  execute(domains, processes);
  domains.erase(domains.begin()+4); // remove level set associated with dummy gate

  // remove SiGe interlayer
  {
    processes.clear();
    auto SiGeEtch = lsSmartPointer<isotropic>::New(std::vector<double>({-0.65, -0.05, -1.0, -0.05}));
    processes.push_back(Process("SiGe-Removal", 10, SiGeEtch));
  }
  execute(domains, processes);
  domains.erase(domains.begin()+2); // Erase level set associated with SiGe

  // gate dielectric, gate metal, gate contact deposition
  {
    processes.clear();
    processes.push_back(Process("HfO2-Deposition", 2, isoVelocity, true));
    processes.push_back(Process("TiN-Deposition", 4, isoVelocity, true));
    processes.push_back(Process("PolySi-Deposition", 20, isoVelocity, true));
  }
  execute(domains, processes);

  // gate CMP
  planarise(domains, 47.5);

  return 0;
}
