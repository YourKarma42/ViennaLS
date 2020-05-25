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


#include <unordered_map>
//#include <hrleDomain.hpp>




//#include <hrleSparseBoxIterator.hpp>
//#include <hrleVectorType.hpp>

#include<lsCalculateCurvatures.hpp>



//____________testing not necessary_________________

#include <chrono>
#include <lsCalculateNormalVectors.hpp>

//____________testing end___________________________

constexpr int D = 2;
typedef double NumericType;


void createPointCloudOutput(lsDomain<NumericType, D> &passedlsDomain, std::string outputName){

  std::cout << "Calculating normal vectors..." << std::endl;




}






lsDomain<double, D> makeSphere(double gridDelta){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0., 0., 0.};
    double radius = 10.0;

    lsDomain<double,D> levelSet(gridDelta);

    lsMakeGeometry<double, D>(levelSet, lsSphere<double, D>(origin, radius)).apply();

      lsMesh narrowband;
      std::cout << "Extracting narrowband..." << std::endl;
      lsToSurfaceMesh<NumericType, D>(levelSet, narrowband).apply();

      lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "directly_after").apply();



    return levelSet;

}


  NumericType calcDist(hrleSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>& starStencil, bool inside){

    //TODO: write eikonal equation somwhere

    //debug
    NumericType tmp_debug[6];

    NumericType stencilMin[D];

    int numUndefined = 0;




    //find the maximum values in the stencil for the Eikonal equation
    for(int i = 0; i < D; i++){

      //TODO: remove v
      NumericType v;

      NumericType pos = starStencil.getNeighbor(i).getValue();

      NumericType neg = starStencil.getNeighbor(i+D).getValue();

      //change the sign to make fast marching on the inside correct
      if(inside){

        pos = pos * -1.;

        neg = neg * -1.;

      }


      tmp_debug[i] = pos;
      tmp_debug[i+3] = neg;

      //check if the current point in the stencil is defined (not +/- inf)
      if( pos !=  lsDomain<NumericType, D>::POS_VALUE && 
          pos !=  lsDomain<NumericType, D>::NEG_VALUE){

        if(neg !=  lsDomain<NumericType, D>::POS_VALUE && 
           neg !=  lsDomain<NumericType, D>::NEG_VALUE){
            if(pos < neg){
              v = pos;
            }else{
              v = neg;
            }
        }else{
          v = pos;
        }

      }else if(neg !=  lsDomain<NumericType, D>::POS_VALUE && 
               neg !=  lsDomain<NumericType, D>::NEG_VALUE){
        v = neg;
      }else{
        numUndefined++;
        v = 0.;
      }

      stencilMin[i] = v; 

    }

    //all points in the stencil are undefined return 0.
    if(numUndefined == D)
      return 0.;

    //if(inside)
    //  return -1.;


    NumericType sol = 0.;

  //debug
 //   for(int i = 0; i< D;i++){
 //    if(stencilMax[i]==0.)
 //       std::cout << "blbu";
 //   }
    
    for(int run = D; run > 0; run--){
      
      NumericType b = 0.;
      NumericType c = (NumericType)run;

      //calculate upwind difference
      for(int i = 0; i < D; i++){
        b += stencilMin[i];
        c -= (NumericType)run * stencilMin[i] * stencilMin[i]; 
      }

      NumericType discriminant = (b*b) + c;

      //discriminant is negative must perform a lower dimensional update
      if(discriminant < 0.){
        //remove max value
        NumericType max = 0.;
        int maxIndex = 0;
        for(int i = 0; i < D; i++){
          if(stencilMin[i] > max){
            max = stencilMin[i];
            maxIndex = i;
          }
        }
        stencilMin[maxIndex] = 0.;
      }else{

        //we are interested in the greater solution of the equation - dose not need to be considered
        sol = (b + std::sqrt(discriminant))/(NumericType)run;
        break;
      }
    }
    

    if(inside){
      return -sol;
    }else{
      return sol;
    }
    

  }


  void testFMM(lsDomain<NumericType, D> * levelSet,
              std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints){

          NumericType incr = 0.;
    //TODO: DO checks of the levelset

    //Do FMM until the narrowband converges
    //The convergence of the narrowband depends on the width that is required
    //TODO: at the moment random static number think of better war to end the loop
    for(int runs=0; runs < 50; runs++){
      
      //debug
      std::vector<NumericType> prevVal;
      std::vector<NumericType> error;
      auto &grid = levelSet->getGrid();

      lsDomain<NumericType, D> newlsDomain(grid);
      typename lsDomain<NumericType, D>::DomainType &newDomain = newlsDomain.getDomain();
      typename lsDomain<NumericType, D>::DomainType &domain = levelSet->getDomain();

      std::cout << "num points: " << domain.getNumberOfPoints() << std::endl;

      //TODO: think of correct value not 2
      //TODO: eventuell leer callen
      //newDomain.initialize(domain.getNewSegmentation(),
      //                      domain.getAllocation() * 4);
      newDomain.initialize();

      int p = 0;

      auto &domainSegment = newDomain.getDomainSegment(p);

      hrleVectorType<hrleIndexType, D> startVector =
          (p == 0) ? grid.getMinGridPoint()
              : newDomain.getSegmentation()[p - 1];

      hrleVectorType<hrleIndexType, D> endVector =
          (p != static_cast<int>(newDomain.getNumberOfSegments() - 1))
              ? newDomain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint()); 

      //TODO: debug
      int count = 0;

      //Fast Marching
      for (hrleSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
        neighborIt(domain, startVector);
        neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &centerIt = neighborIt.getCenter();


        //mark all lvl set values as accepted
        if((lsPoints.find(centerIt.getStartIndices()) != lsPoints.end())){
          //add them to the new lvl set
          domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                              centerIt.getValue()); 
          //debug
          NumericType radius = 0.;
          
          for(int i = 0; i < D; i++){
            NumericType coord = centerIt.getStartIndices(i) * grid.getGridDelta();
            radius += coord * coord;
          }
          error.push_back(std::abs(10. - std::sqrt(radius) - centerIt.getValue()));

          prevVal.push_back(0.);   
          //debug                                         
          count++;

        }else{    

          NumericType old = centerIt.getValue();

         /* if(std::abs(old) <= 0.5 && runs >= 8){
            std::cout << "blub" << std::endl;

          }*/

          //calculate distance to the interface using the iconal equation 
          NumericType dist = calcDist(neighborIt, (centerIt.getValue() < 0.));
       
          //check if a distance could be calculated for the current point
          if(dist != 0.){
            domainSegment.insertNextDefinedPoint(neighborIt.getIndices(), 
            //(centerIt.getValue()<0.) ? -dist : dist);
            dist);

            //debug
            if((old != lsDomain<NumericType, D>::POS_VALUE 
             && old != lsDomain<NumericType, D>::NEG_VALUE)){
              prevVal.push_back(std::abs(dist-old));

              NumericType radius = 0.;

              for(int i = 0; i < D; i++){
                NumericType coord = centerIt.getStartIndices(i) * grid.getGridDelta();
                radius += coord * coord;
              }
              error.push_back(std::abs(10. - std::sqrt(radius)- dist));
            }else{
              prevVal.push_back(0.);
              error.push_back(10.);
            }
          }else{
            domainSegment.insertNextUndefinedPoint(neighborIt.getIndices(), 
            (centerIt.getValue()<0.) ? lsDomain<NumericType, D>::NEG_VALUE : lsDomain<NumericType, D>::POS_VALUE);

          }
          
        }
      }

      newDomain.finalize();
      levelSet->deepCopy(newlsDomain);


      //debug

      if(D<3){
        //output takes to long in 3d
        lsMesh narrowband;
        std::cout << "Extracting narrowband..." << std::endl;
        lsToMesh<NumericType, D>(*levelSet, narrowband, true, false).apply();

        narrowband.insertNextScalarData(prevVal, "difference");

        narrowband.insertNextScalarData(error, "error");

        lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "FM_norm_S" + std::to_string(runs) ).apply();

        std::cout << "kept points: " << count << std::endl;
      }




    } 

    //TODO: constant is stupid change to dynamic value e.g stop condition or fixed amount of steps
    NumericType width = 2.;
    levelSet->getDomain().segment();
    levelSet->finalize(width);
      
  
  }

  lsDomain<NumericType, D> testConvert(lsDomain<NumericType, D> &passedlsDomain){


    std::cout << "Converting Distances..." << std::endl;

    //increase lvl set size
    int order = 1;

    //test!

    //lsExpand<NumericType, D>(passedlsDomain, 2 * (order + 2) + 1).apply();


    double pointsPerSegment =
    double(2 * passedlsDomain.getDomain().getNumberOfPoints()) /
    double(passedlsDomain.getLevelSetWidth());

    auto grid = passedlsDomain.getGrid();

    //! Calculate Normalvectors
    int p = 0;
#ifdef _OPENMP
    p = omp_get_thread_num();
#endif

    hrleVectorType<hrleIndexType, D> startVector =
        (p == 0) ? grid.getMinGridPoint()
                : passedlsDomain.getDomain().getSegmentation()[p - 1];

    hrleVectorType<hrleIndexType, D> endVector =
      (p != static_cast<int>(passedlsDomain.getNumberOfSegments() - 1))
        ? passedlsDomain.getDomain().getSegmentation()[p]
        : grid.incrementIndices(grid.getMaxGridPoint());

    //calculation of normal
    const NumericType gridDelta = passedlsDomain.getGrid().getGridDelta();

    typedef typename lsDomain<NumericType, D>::PointValueVectorType pointDataType;

    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints;
    
    lsPoints.reserve(passedlsDomain.getNumberOfPoints());
    

    hrleVectorType<hrleIndexType, D> index;


    //find aktive ls points
    for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
            neighborIt(passedlsDomain.getDomain(), startVector);
           neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &center = neighborIt.getCenter();
        if (!center.isDefined()) {
          continue;
        } else if (std::abs(center.getValue()) > 0.5) {
          // undefined run
          continue;
        }
      
        lsPoints.insert(center.getStartIndices());
    }

    testFMM(& passedlsDomain , lsPoints);

    //lsExpand<NumericType, D>(passedlsDomain, 2 * (order + 2) + 1).apply();

    return passedlsDomain;

        
}

  lsDomain<NumericType, D> ConvertLS(lsDomain<NumericType, D> &passedlsDomain){


    std::cout << "Converting Distances..." << std::endl;

    //increase lvl set size
    int order = 1;

    lsExpand<NumericType, D>(passedlsDomain, 2 * (order + 2) + 1).apply();


    double pointsPerSegment =
    double(2 * passedlsDomain.getDomain().getNumberOfPoints()) /
    double(passedlsDomain.getLevelSetWidth());

    auto grid = passedlsDomain.getGrid();

    typename lsDomain<NumericType, D>::DomainType &oldDomain = passedlsDomain.getDomain();

    const NumericType gridDelta = passedlsDomain.getGrid().getGridDelta();

    //Initialize new level set with normalized values
    lsDomain<NumericType,D> newLS(gridDelta);

    auto newGrid = newLS.getGrid();

    typename lsDomain<NumericType, D>::DomainType &newDomain = newLS.getDomain();

    newDomain.initialize();


    //! Calculate Normalvectors
    int p = 0;
#ifdef _OPENMP
    p = omp_get_thread_num();
#endif

    auto &newDomainSegment = newDomain.getDomainSegment(p);

    hrleVectorType<hrleIndexType, D> startVector =
        (p == 0) ? grid.getMinGridPoint()
                : passedlsDomain.getDomain().getSegmentation()[p - 1];

    hrleVectorType<hrleIndexType, D> endVector =
      (p != static_cast<int>(passedlsDomain.getNumberOfSegments() - 1))
        ? passedlsDomain.getDomain().getSegmentation()[p]
        : grid.incrementIndices(grid.getMaxGridPoint());


    //TODO: probalby remove dont needed keep it for now
    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints;
    
    std::vector<double> oldVal;

    oldVal.reserve(passedlsDomain.getNumberOfPoints());
    lsPoints.reserve(passedlsDomain.getNumberOfPoints());
    

    hrleVectorType<hrleIndexType, D> index;

    //debug
    lsMesh test_output;
    test_output.clear();
    std::vector<double> lsVals;

    lsMesh test_output2;
    test_output2.clear();
    std::vector<double> lsVals2;


    for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
          neighborIt(passedlsDomain.getDomain(), startVector);
          neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &centerIt = neighborIt.getCenter();
        if (!centerIt.isDefined() || std::abs(centerIt.getValue()) > 0.5) {
          //write undefined run in new level set
          newDomain.getDomainSegment(p).insertNextUndefinedPoint(neighborIt.getIndices(), 
          (centerIt.getValue()<0.) ? lsDomain<NumericType, D>::NEG_VALUE : lsDomain<NumericType, D>::POS_VALUE);

          continue;
        } 

        
        lsPoints.insert(centerIt.getStartIndices());


        //TODO: debug save old LS value for comparison later
        NumericType oldLSValue = centerIt.getValue();

        //calculate normal vector of defined grid point
        std::array<NumericType, D> n;
        NumericType normN = 0.;
        for (int i = 0; i < D; i++) {
          
          if(centerIt.getValue() == 0.){
       //     std::cout << "p " << neighborIt.getNeighbor(i).getValue() << std::endl;
       //     std::cout << "n " << neighborIt.getNeighbor(i + D).getValue() << std::endl;

          }

          NumericType pos1 = neighborIt.getNeighbor(i).getValue();
          NumericType neg1 = neighborIt.getNeighbor(i + D).getValue();
          n[i] = (pos1 - neg1) * 0.5;
          
          normN += n[i] * n[i];
        }
        normN = std::sqrt(normN);

        for (int i = 0; i < D; i++) {
          n[i] /= normN;
        }

        //renormalized ls value in direction of the normal vector
        NumericType max = 0.;
        for (unsigned i = 0; i < D; ++i) {

            // grid position
            index[i] = centerIt.getStartIndices(i);

            if (std::abs(n[i]) > max) {
            max = std::abs(n[i]);
            }
        }

        //development
        NumericType normSurfPoint = 0.;
        NumericType normGridPoint = 0.;

        std::array<double, 3> node;
        node[2] = 0.;
        NumericType scaling = centerIt.getValue() * gridDelta * max;
        for (unsigned i = 0; i < D; ++i) {
          node[i] = double(centerIt.getStartIndices(i)) * gridDelta;

          normGridPoint += node[i] * node[i];

          //rescale vector to touch the surface
          node[i] -= scaling * n[i];

          normSurfPoint += node[i] * node[i];
        }

        normSurfPoint = std::sqrt(normSurfPoint);
        normGridPoint = std::sqrt(normGridPoint);

  //      if(centerIt.getValue()<0.0)
 //         std::cout << "blbu";

        std::array<unsigned, 1> vertex;
        vertex[0] =test_output.nodes.size();
        test_output.insertNextVertex(vertex);

       // node[2] = 0.;
       // for (unsigned i = 0; i < D; ++i) {
       //   node[i] -= scaling * n[i];
       // }

        test_output.insertNextNode(node);

        lsVals.push_back(scaling);

        NumericType test = normGridPoint-normSurfPoint;

        //std::cout << "subtr:" << test << std::endl;
        //std::cout << "scale:" << scaling << std::endl;

        vertex[0] =test_output2.nodes.size();
        test_output2.insertNextVertex(vertex);

        node[2] = 0.;
        for (unsigned i = 0; i < D; ++i) {
          node[i] = double(centerIt.getStartIndices(i)) * gridDelta;
        }

        test_output2.insertNextNode(node);

        lsVals2.push_back(oldLSValue);


              
        NumericType newLsValue = scaling;
        //NumericType newLsValue = ((centerIt.getValue() < 0.) ? -1. : 1.) * test;



        newDomain.getDomainSegment(p).insertNextDefinedPoint(neighborIt.getIndices(), newLsValue);

        oldVal.push_back(centerIt.getValue());

    }

    newDomain.finalize();

    test_output.insertNextScalarData(lsVals, "new Vals");
    test_output.insertNextScalarData(lsVals2, "old Vals");

    lsVTKWriter(test_output, lsFileFormatEnum::VTU , "gridNew" ).apply();

    test_output2.insertNextScalarData(lsVals, "new Vals"); 
    test_output2.insertNextScalarData(lsVals2, "old Vals");

    lsVTKWriter(test_output2, lsFileFormatEnum::VTU , "gridOld" ).apply();


//  newLS  passedlsDomain
    testFMM(& newLS , lsPoints);

    
    return newLS;

   

       testFMM(& passedlsDomain , lsPoints);

    //debug
    //return newLS;
    return passedlsDomain;

        
}


  lsDomain<NumericType, D> dontConvert(lsDomain<NumericType, D> &passedlsDomain){

    std::cout << "NOT Converting Distances..." << std::endl;

    auto grid = passedlsDomain.getGrid();

    //! Calculate Normalvectors
    int p = 0;
#ifdef _OPENMP
    p = omp_get_thread_num();
#endif

    hrleVectorType<hrleIndexType, D> startVector =
        (p == 0) ? grid.getMinGridPoint()
                : passedlsDomain.getDomain().getSegmentation()[p - 1];

    hrleVectorType<hrleIndexType, D> endVector =
      (p != static_cast<int>(passedlsDomain.getNumberOfSegments() - 1))
        ? passedlsDomain.getDomain().getSegmentation()[p]
        : grid.incrementIndices(grid.getMaxGridPoint());

    const NumericType gridDelta = passedlsDomain.getGrid().getGridDelta();


    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints;
    

    for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
            neighborIt(passedlsDomain.getDomain(), startVector);
           neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &center = neighborIt.getCenter();
        if (!center.isDefined()) {
          continue;
        } else if (std::abs(center.getValue()) > 0.5) {
          
          // undefined run
          continue;
        }

        
        lsPoints.insert(center.getStartIndices());

    }


//  newLS  passedlsDomain
    testFMM(& passedlsDomain , lsPoints);

      lsMesh narrowband;
      std::cout << "Extracting narrowband..." << std::endl;
      lsToMesh<NumericType, D>(passedlsDomain, narrowband, true, false).apply();

      lsVTKWriter(narrowband, lsFileFormatEnum::VTU ,  "AfterNotConverting").apply();


    return passedlsDomain;
   

        
}




int main() {


    omp_set_num_threads(1);

    NumericType gridDelta = 0.25;


    std::vector<lsDomain<NumericType, D> *> levelSets;

    lsDomain<NumericType,D> levelSet = makeSphere(gridDelta);

    lsMesh narrowband1;

    lsToMesh<NumericType, D>(levelSet, narrowband1, true, true).apply();

    lsVTKWriter(narrowband1, lsFileFormatEnum::VTU , "OriginalSphere" ).apply();

    lsMesh diskmesh;

    lsToDiskMesh<NumericType, D>(levelSet, diskmesh).apply();

    lsVTKWriter(diskmesh, lsFileFormatEnum::VTU , "SphereDiskMesh" ).apply();

    levelSets.push_back(&levelSet);  

    //Converts the LS from sparse manhatten normalization to Euklid normalization
    //lsDomain<NumericType,D> newLevelSet =  ConvertLS(*(levelSets.back()));

    lsDomain<NumericType,D> newLevelSet =  dontConvert(*(levelSets.back()));

   

    //std::cout << newLevelSet.getNumberOfSegments() << std::endl;

    //passes on the created sphere
    //lsDomain<NumericType,D> newLevelSet =  testConvert(*(levelSets.back()));

    //perform FMM

    //calculate curvatures

    //write output

    //lsCalculateNormalVectors<NumericType, D> test_normals(newLevelSet);

    //test_normals.apply();

    //auto& normalVectors = newLevelSet.getNormalVectors();

    std::vector<std::array<NumericType, 3>> normal;

    std::vector<NumericType> grad;


   
   //     ? passedlsDomain.getDomain().getSegmentation()[p]
   //     : grid.incrementIndices(grid.getMaxGridPoint());

    for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
          neighborIt(newLevelSet.getDomain(), newLevelSet.getGrid().getMinGridPoint());
          neighborIt.getIndices() < newLevelSet.getGrid().incrementIndices(newLevelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

        auto &centerIt = neighborIt.getCenter();
        if (!centerIt.isDefined() || std::abs(centerIt.getValue()) > 0.5) {
          continue;
        } 

        //calculate normal vector of defined grid point
        std::array<NumericType, 3> n;
        n[2] = 0.;
        for (int i = 0; i < D; i++) {
          NumericType pos1 = neighborIt.getNeighbor(i).getValue();
          NumericType neg1 = neighborIt.getNeighbor(i + D).getValue();
          n[i] = (pos1 - neg1) /(2.);
          
        }

        NumericType normGrad = 0.;

        for(int i = 0; i < D; i++){
          normGrad += n[i]*n[i];
        }

        grad.push_back(std::sqrt(normGrad));
              
        normal.push_back(n);



    }

    lsMesh narrowband;
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(newLevelSet, narrowband, true, true).apply();


    narrowband.insertNextVectorData(normal, "Normals");
    narrowband.insertNextScalarData(grad, "norm of grad");

    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "final_output" ).apply();

    lsMesh mesh;
    std::cout << "Extracting surface mesh..." << std::endl;
    lsToSurfaceMesh<double, D>(newLevelSet, mesh).apply();


    lsVTKWriter(mesh, lsFileFormatEnum::VTU ,"final_output_mesh").apply();

    std::cout << "Finished" << std::endl;
 

    return 0;
}



/*
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
*/

