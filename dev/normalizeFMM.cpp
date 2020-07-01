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


#include <unordered_set>
//#include <hrleDomain.hpp>




//#include <hrleSparseBoxIterator.hpp>
//#include <hrleVectorType.hpp>

#include <lsCalculateCurvatures.hpp>

#include <omp.h>

#include <../dev/derivatives.hpp>



//____________testing not necessary_________________

#include <chrono>
#include <lsCalculateNormalVectors.hpp>

//____________testing end___________________________

constexpr int D = 3;
typedef double NumericType;



lsDomain<double, D> makeSphere(double gridDelta){

    std::cout << "creating sphere..." << std::endl;

    double origin[3] = {0., 0., 0.};
    double radius = 10.0;

    lsDomain<double,D> levelSet(gridDelta);

    lsMakeGeometry<double, D>(levelSet, lsSphere<double, D>(origin, radius)).apply();


    return levelSet;

}


  NumericType calcDist(hrleSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>& starStencil, bool inside){

    //TODO: write eikonal equation somwhere


    NumericType stencilMin[D];

    int numUndefined = 0;


    //find the maximum values in the stencil for the Eikonal equation
    for(int i = 0; i < D; i++){


      NumericType pos = starStencil.getNeighbor(i).getValue();

      NumericType neg = starStencil.getNeighbor(i+D).getValue();

      //change the sign to make fast marching on the inside correct
      if(inside){
        pos = pos * -1.;
        neg = neg * -1.;
      }

      //check if the current point in the stencil is defined (not +/- inf)
      if( pos !=  lsDomain<NumericType, D>::POS_VALUE && 
          pos !=  lsDomain<NumericType, D>::NEG_VALUE){

        if(neg !=  lsDomain<NumericType, D>::POS_VALUE && 
           neg !=  lsDomain<NumericType, D>::NEG_VALUE){
            if(pos < neg){
              stencilMin[i] = pos;
            }else{
              stencilMin[i] = neg;
            }
        }else{
          stencilMin[i] = pos;
        }

      }else if(neg !=  lsDomain<NumericType, D>::POS_VALUE && 
               neg !=  lsDomain<NumericType, D>::NEG_VALUE){
        stencilMin[i] = neg;
      }else{
        numUndefined++;
        stencilMin[i] = 0.;
      }

    }

    //all points in the stencil are undefined return 0.
    if(numUndefined == D)
      return 0.;

    NumericType sol = 0.;
    
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
    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & lsPoints){
              //,){

    std::cout << "Doing FM..." << std::endl;


    NumericType incr = 0.;
    //TODO: DO checks of the levelset

    //Do FMM until the narrowband converges
    //The convergence of the narrowband depends on the width that is required
    //TODO: at the moment random static number think of better war to end the loop
    for(int runs=0; runs < 20; runs++){
      
      //debug
      std::vector<NumericType> prevVal;
      std::vector<NumericType> error;
      std::vector<NumericType> radiusPoint;
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

          radius = std::sqrt(radius);

          radiusPoint.push_back(radius);

          if(centerIt.getValue() >= 0.){
              error.push_back( std::abs((std::abs(centerIt.getValue())*grid.getGridDelta()) - radius) );
          }else{
              error.push_back( std::abs((std::abs(centerIt.getValue())*grid.getGridDelta()) + radius) );
          }

          prevVal.push_back(0.);   
          //debug                                         
          count++;

        }else{    

          NumericType old = centerIt.getValue();

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

              radius = std::sqrt(radius);

              radiusPoint.push_back(radius);

              if(radius >= 10.){
                error.push_back( std::abs((std::abs(dist)*grid.getGridDelta()) - radius) );


              }else{
                error.push_back( std::abs((std::abs(dist)*grid.getGridDelta()) + radius) );

              }


            }else{
              prevVal.push_back(0.);

              if(old == lsDomain<NumericType, D>::POS_VALUE){
                error.push_back(100.);
              }else{
                error.push_back( 100.);
              }

              radiusPoint.push_back(100.);
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

      if(D<3 || runs==19){
        //output takes to long in 3d
        lsMesh narrowband;
        std::cout << "Extracting narrowband..." << std::endl;
        lsToMesh<NumericType, D>(*levelSet, narrowband, true, false).apply();

        narrowband.insertNextScalarData(prevVal, "difference");

        narrowband.insertNextScalarData(error, "error");

        narrowband.insertNextScalarData(radiusPoint, "Radius");

        lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "FM_norm_S" + std::to_string(runs) ).apply();

        std::cout << "kept points: " << count << std::endl;
      }




    } 

    //TODO: constant is stupid change to dynamic value e.g stop condition or fixed amount of steps
    NumericType width = 3.;
    levelSet->getDomain().segment();
    levelSet->finalize(width);
      
  
  }


  lsDomain<NumericType, D> ConvertLS(lsDomain<NumericType, D> &passedlsDomain, 
    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & lsPoints){

    //TODO:remove
    std::cout << "Converting Distances..." << std::endl;

    //increase lvl set size
    lsExpand<NumericType, D>(passedlsDomain, 7).apply();


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

    newDomain.initialize(oldDomain.getNewSegmentation(), oldDomain.getAllocation());


    //! Calculate Normalvectors
    int p = 0;
#ifdef _OPENMP
    p = omp_get_thread_num();
#endif

    auto &newDomainSegment = newDomain.getDomainSegment(p);

    hrleVectorType<hrleIndexType, D> startVector =
        (p == 0) ? grid.getMinGridPoint()
                : newDomain.getSegmentation()[p - 1];

    hrleVectorType<hrleIndexType, D> endVector =
      (p != static_cast<int>(newDomain.getNumberOfSegments() - 1))
        ? newDomain.getSegmentation()[p]
        : grid.incrementIndices(grid.getMaxGridPoint());




    //TODO: remove
    std::vector<double> oldVal;

    oldVal.reserve(passedlsDomain.getNumberOfPoints());
    //lsPoints.reserve(passedlsDomain.getNumberOfPoints());
    

    hrleVectorType<hrleIndexType, D> index;

 
    for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
          neighborIt(passedlsDomain.getDomain(), startVector);
          neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &centerIt = neighborIt.getCenter();
        if (!centerIt.isDefined() || (std::abs(centerIt.getValue()) > 0.5) ) {
          //write undefined run in new level set
          newDomain.getDomainSegment(p).insertNextUndefinedPoint(neighborIt.getIndices(), 
          (centerIt.getValue()<0.) ? lsDomain<NumericType, D>::NEG_VALUE : lsDomain<NumericType, D>::POS_VALUE);

          continue;
        } 

        //TODO: probably remove
        lsPoints.insert(centerIt.getStartIndices());


        //TODO: debug save old LS value for comparison later
        NumericType oldLSValue = centerIt.getValue();

        //calculate normal vector of defined grid point
        std::array<NumericType, D> n;
        NumericType normN = 0.;
        for (int i = 0; i < D; i++) {
          
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
   
        NumericType newLsValue = centerIt.getValue() * max;

        newDomain.getDomainSegment(p).insertNextDefinedPoint(neighborIt.getIndices(), newLsValue);

        oldVal.push_back(centerIt.getValue());

    }

    newDomain.finalize();


    return newLS;
        
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


    return passedlsDomain;
        
}

  lsDomain<NumericType, D> hardConvertCircle(lsDomain<NumericType, D> &passedlsDomain,
    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & lsPoints){

    std::cout << "Hard converting Circle..." << std::endl;

    auto grid = passedlsDomain.getGrid();

    NumericType gridDelta = grid.getGridDelta();

    lsDomain<NumericType,D> newLS(gridDelta);

    auto newGrid = newLS.getGrid();

    typename lsDomain<NumericType, D>::DomainType &newDomain = newLS.getDomain();

    newDomain.initialize(passedlsDomain.getDomain().getNewSegmentation(), passedlsDomain.getDomain().getAllocation());

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
   

    for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
            neighborIt(passedlsDomain.getDomain(), startVector);
           neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &center = neighborIt.getCenter();

        
        if (!center.isDefined() ||  std::abs(center.getValue()) > 0.5) {

          newDomain.getDomainSegment(p).insertNextUndefinedPoint(neighborIt.getIndices(), 
          (center.getValue()<0.) ? lsDomain<NumericType, D>::NEG_VALUE : lsDomain<NumericType, D>::POS_VALUE);
          
          // undefined run
          continue;
        }
        NumericType radius = 0.;
        for(int i = 0; i < D; i++){
          NumericType coord = center.getStartIndices(i) * gridDelta;
          radius += coord * coord;
        }

        if(std::sqrt(radius) >= 10){
          newLS.getDomain().getDomainSegment(p).insertNextDefinedPoint(neighborIt.getIndices(), std::sqrt(radius) - 10.);
        }else{
          newLS.getDomain().getDomainSegment(p).insertNextDefinedPoint(neighborIt.getIndices(), 10. - std::sqrt(radius));
        }
        
        lsPoints.insert(center.getStartIndices());

    }


    return newLS;
        
}

NumericType convertToFuncVal(NumericType lsVal, hrleVectorType<hrleIndexType, D> gridPoint, NumericType gridDelta){

  NumericType distGridPoint = 0.;

  for(int i = 0; i < D; i++){
    distGridPoint += (gridPoint[i]*gridDelta) * (gridPoint[i]*gridDelta);
  }

  distGridPoint = std::sqrt(distGridPoint);

  if(lsVal > 0.){
    return (lsVal*gridDelta) + distGridPoint;
  }else{
    return ((lsVal*gridDelta) + distGridPoint);
  }



}

NumericType calc_mean_curvature(std::vector<NumericType> d){


    //TODO: write used formula and paper with reference
    //The 2 in the denominator is important when considering the mean curvature in the context of differential geometry in books for
    //lvl-set functions the 1/2 is often missing


    NumericType norm_grad_pow3 = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;

    //expanded fom of the equation
    return
    //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx + F_yy)
    (d[0]*d[0]*(d[4] + d[5]) + d[1]*d[1]*(d[3] + d[5]) + d[2]*d[2]*(d[3] + d[4]) +

    //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
    -2.*(d[0]*d[1]*d[6] + d[0]*d[2]*d[8] + d[1]*d[2]*d[7]))
            
    // /2*(F_x² + F_y² + F_z²)^(3/2)
    /(2.*norm_grad_pow3);

      //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)

}


void create_output(lsDomain<NumericType,D> & levelSet,
  std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & lsPoints){

    hrleSparseBoxIterator<hrleDomain<NumericType, D>> neighborIterator(levelSet.getDomain(), 2);

     std::vector<std::array<NumericType, 3>> normal;

    std::vector<std::array<NumericType, 3>> secondOrderDerivatives1;

    std::vector<std::array<NumericType, 3>> secondOrderDerivatives2;



    std::vector<NumericType> grad1;

    std::vector<NumericType> meanCurvature1;

    std::vector<NumericType> meanCurvatureSD1;


    std::vector<NumericType> grad2;

    std::vector<NumericType> meanCurvature2;

    std::vector<NumericType> meanCurvatureSD2;

    std::vector<NumericType> meanCurvature3;

    std::vector<NumericType> meanCurvatureNew;

    NumericType gridDelta = levelSet.getGrid().getGridDelta();

    curvatur1<NumericType, D> calculatorTest(gridDelta);







   
   //     ? passedlsDomain.getDomain().getSegmentation()[p]
   //     : grid.incrementIndices(grid.getMaxGridPoint());

    for (hrleConstSparseStarIterator<typename lsDomain<NumericType, D>::DomainType>
          neighborIt(levelSet.getDomain(), levelSet.getGrid().getMinGridPoint());
          neighborIt.getIndices() < levelSet.getGrid().incrementIndices(levelSet.getGrid().getMaxGridPoint()); neighborIt.next()) {

        auto &centerIt = neighborIt.getCenter();
        if (!centerIt.isDefined() || (lsPoints.find(centerIt.getStartIndices()) == lsPoints.end())) {
          continue;
        } 

        // move neighborIterator to current position
        neighborIterator.goToIndicesSequential(centerIt.getStartIndices());

        



        //NumericType gridDelta = levelSet.getGrid().getGridDelta();

        //calculate normal vector of defined grid point
        std::array<NumericType, 3> n;
        n[2] = 0.;
        std::array<NumericType, 3> d1;
        d1[2] = 0.;
        std::array<NumericType, 3> d2;
        d2[2] = 0.;

        //std::cout << centerIt.getStartIndices() << std::endl;

        std::vector<NumericType> derivatives1(9);

        std::vector<NumericType> derivatives2(9);

        std::vector<NumericType> derivatives3(9);



        std::vector<NumericType> d_p(3);

        std::vector<NumericType> d_n(3);

        std::vector<NumericType> d_c(3);

        std::vector<NumericType> a;



        for (int i = 0; i < D; i++) {

          hrleVectorType<hrleIndexType, D> posUnit(0);
          hrleVectorType<hrleIndexType, D> negUnit(0);

          hrleVectorType<hrleIndexType, D> test = neighborIterator.getIndices();

          posUnit[i] = 1;
          negUnit[i] = -1;

          int second_pos = (i+1) % D;

          if(i == 1){
            second_pos = 0;
          }

          //get required ls values
          NumericType phi_0 = neighborIterator.getCenter().getValue();

          phi_0 = convertToFuncVal(phi_0, neighborIterator.getCenter().getStartIndices(), gridDelta);

          NumericType phi_px = neighborIterator.getNeighbor(posUnit).getValue();
          NumericType phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

          phi_px = convertToFuncVal(phi_px, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
          phi_nx = convertToFuncVal(phi_nx, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

          posUnit[i] = 2;
          negUnit[i] = -2;

          NumericType phi_p2x = neighborIterator.getNeighbor(posUnit).getValue();
          NumericType phi_n2x = neighborIterator.getNeighbor(negUnit).getValue();

          phi_p2x = convertToFuncVal(phi_p2x, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
          phi_n2x = convertToFuncVal(phi_n2x, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

          posUnit[i] = 1;
          negUnit[i] = -1;

          posUnit[second_pos] = 1;
          negUnit[second_pos] = 1;

          NumericType phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
          NumericType phi_np = neighborIterator.getNeighbor(negUnit).getValue();

          phi_pp = convertToFuncVal(phi_pp, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
          phi_np = convertToFuncVal(phi_np, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

          posUnit[second_pos] = -1;
          negUnit[second_pos] = -1;

          NumericType phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
          NumericType phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

          phi_pn = convertToFuncVal(phi_pn, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
          phi_nn = convertToFuncVal(phi_nn, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

          posUnit[i] = 0;
          negUnit[i] = 0;

          posUnit[second_pos] = 1;
          negUnit[second_pos] = -1;

          NumericType phi_py = neighborIterator.getNeighbor(posUnit).getValue();
          NumericType phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

          phi_py = convertToFuncVal(phi_py, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
          phi_ny = convertToFuncVal(phi_ny, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

          //test

          /*d_p[i] = phi_px - phi_0;

          d_n[i] = phi_0 - phi_nx;

          d_c[i] = (phi_px - phi_nx)*0.5;

          phi_pp - phi_np;

          phi_pn - phi_nn;*/





          // 1 1 1
          derivatives1[i] = (phi_px - phi_nx)/(2.*gridDelta);

          derivatives1[i+3] = (phi_px - 2.*phi_0 + phi_nx)/(gridDelta*gridDelta);

          derivatives1[i+6] = (phi_pp - phi_pn -phi_np + phi_nn)/(4.*gridDelta*gridDelta);


          // 2 2 2 
          derivatives2[i] = (phi_pp - phi_np + phi_pn - phi_nn)/(4.*gridDelta);

          derivatives2[i+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)/(gridDelta*gridDelta*3.0);

          derivatives2[i+6] = -(phi_px + phi_nx + phi_py + phi_ny - 2. * phi_0 - phi_pp - phi_nn)/(2.*gridDelta*gridDelta);
          


          // 2 2 1
          derivatives3[i] = (phi_pp - phi_np + phi_pn - phi_nn)/(4.*gridDelta);

          derivatives3[i+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)/(gridDelta*gridDelta*3.0);

          derivatives3[i+6] = (phi_pp - phi_pn -phi_np + phi_nn)/(4.*gridDelta*gridDelta);


          //second order F_xx
          //d2[i] = (-phi_p2x + 16.*phi_px - 30.*phi_0 + 16.*phi_nx - phi_n2x)/(12.*gridDelta*gridDelta);
          
        }


        NumericType normGrad1 = 0.;

        NumericType mCurve1 = 0.;

        NumericType normGrad2 = 0.;

        NumericType mCurve2 = 0.;

        for(int i = 0; i < D; i++){
          normGrad1 += derivatives1[i]*derivatives1[i];
          mCurve1 += derivatives1[i+3];

          normGrad2 += derivatives2[i]*derivatives2[i];
          mCurve2 += derivatives2[i+3];
        }

        meanCurvatureNew.push_back(calculatorTest(neighborIterator));

        grad1.push_back(std::sqrt(normGrad1));

        grad2.push_back(std::sqrt(normGrad2));

        meanCurvatureSD1.push_back((mCurve1 * 0.5));

        meanCurvatureSD2.push_back((mCurve2 * 0.5));

        meanCurvature1.push_back(calc_mean_curvature(derivatives1));

        meanCurvature2.push_back(calc_mean_curvature(derivatives2));

        meanCurvature3.push_back(calc_mean_curvature(derivatives3));
              
        //normal.push_back(n);
        //secondOrderDerivatives1.push_back(d1);
        //secondOrderDerivatives2.push_back(d2);

    }

    lsMesh narrowband;
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(levelSet, narrowband, true, true).apply(lsPoints);
    //lsPoints


    //narrowband.insertNextVectorData(normal, "Normals");
    narrowband.insertNextScalarData(grad1, "norm of grad");
    narrowband.insertNextScalarData(meanCurvatureSD1, "mean c SD");
    narrowband.insertNextScalarData(meanCurvature1, "mean c F");

    //narrowband.insertNextVectorData(normal, "Normals 2");
    narrowband.insertNextScalarData(grad2, "norm of grad 2");
    narrowband.insertNextScalarData(meanCurvatureSD2, "mean c SD 2");
    narrowband.insertNextScalarData(meanCurvature2, "mean c F 2");

    narrowband.insertNextScalarData(meanCurvature3, "mean c F 3");

    narrowband.insertNextScalarData(meanCurvatureNew, "mean new");
    //narrowband.insertNextVectorData(secondOrderDerivatives1, "F_xx_1");
    //narrowband.insertNextVectorData(secondOrderDerivatives2, "F_xx_2");
  
    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "final_output" ).apply();

    lsMesh mesh;
    std::cout << "Extracting surface mesh..." << std::endl;
    lsToSurfaceMesh<double, D>(levelSet, mesh).apply();


    lsVTKWriter(mesh, lsFileFormatEnum::VTU ,"final_output_mesh").apply();

}




int main() {


    omp_set_num_threads(1);

    NumericType gridDelta = 0.25;

    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> lsPoints;


    std::vector<lsDomain<NumericType, D> *> levelSets;

    lsDomain<NumericType,D> levelSet = makeSphere(gridDelta);

    lsMesh diskmesh;

    lsToDiskMesh<NumericType, D>(levelSet, diskmesh).apply();

    lsVTKWriter(diskmesh, lsFileFormatEnum::VTU , "SphereDiskMesh" ).apply();

    levelSets.push_back(&levelSet);  

    test1();

   /* lsExpand(*(levelSets.back()), 2);

    lsCalculateNormalVectors<NumericType, D> test_normals(*(levelSets.back()));

    test_normals.apply();

    lsMesh narrowband1;
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<NumericType, D>(*(levelSets.back()), narrowband1, true, true).apply();

    lsVTKWriter(narrowband1, lsFileFormatEnum::VTU , "standard_normals" ).apply();*/

    //Converts the LS from sparse manhatten normalization to Euklid normalization

    //lsDomain<NumericType,D> newLevelSet =  hardConvertCircle(*(levelSets.back()), lsPoints);

    lsDomain<NumericType,D> newLevelSet =  ConvertLS(*(levelSets.back()), lsPoints);

    //lsDomain<NumericType,D> newLevelSet =  dontConvert(*(levelSets.back()));


    testFMM(& newLevelSet, lsPoints);



    //passes on the created sphere
    //lsDomain<NumericType,D> newLevelSet =  testConvert(*(levelSets.back()));

    create_output(newLevelSet, lsPoints);


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

