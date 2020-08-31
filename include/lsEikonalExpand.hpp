#ifndef LS_EIKONAL_EXPAND_HPP
#define LS_EIKONAL_EXPAND_HPP

//TODO: paralellization

#include <lsPreCompileMacros.hpp>

#include <hrleSparseIterator.hpp>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>

#include <unordered_set>

//TODO: expand!
/// This class expands the level set to a certain numeber of layers by using the Eikonal equation.
/// This expansion is more accurate than using Manhatten expansion, however it has alonger Calculation time.
/// This form of expansion keeps the signed distance porperty of the level set |grad(f)| = 1.

template <class T, int D> class lsEikonalExpand {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  lsDomain<T, D> *levelSet = nullptr;

  T gridDeltaSqared = 0.;

  //TOD: make call by reference (problem with initialization)
  std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash>  activePoints;


public:
  lsEikonalExpand() {}

  lsEikonalExpand(lsDomain<T, D> &passedLevelSet, 
      std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash>  passedActivePoints)
      : levelSet(&passedLevelSet), activePoints(passedActivePoints) {
        //TODO: check if level set is euler normalized and give error if not
      gridDeltaSqared = levelSet->getGrid().getGridDelta() * levelSet->getGrid().getGridDelta();
  }

  void setLevelSet(lsDomain<T, D> &passedLevelSet) {
    levelSet = &passedLevelSet;
  }

  void setActivePoints(std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash>  passedActivePoints) {
    activePoints = passedActivePoints;
  }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsEikonalExpand.")
          .print();
      return;
    }
    //TODO: probably remove
    //T incr = 0.;

    //Do FMM until the narrowband converges
    //The convergence of the narrowband depends on the width that is required
    //TODO: at the moment random static number think of better way to end the loop
    for(int runs=0; runs < 10; runs++){

      const int allocationFactor =
          1 + 1.0 / static_cast<double>(runs);
      //TODO:probably remove
      //const T limit = (runs + 1) * T(0.5);
      
      auto &grid = levelSet->getGrid();
      lsDomain<T, D> newlsDomain(grid);
      typename lsDomain<T, D>::DomainType &newDomain = newlsDomain.getDomain();
      typename lsDomain<T, D>::DomainType &domain = levelSet->getDomain();

      newDomain.initialize(domain.getNewSegmentation(),
                           domain.getAllocation() * allocationFactor);

#pragma omp parallel num_threads(newDomain.getNumberOfSegments())
      {
        int p = 0;
#ifdef _OPENMP
        p = omp_get_thread_num();
#endif

        auto &domainSegment = newDomain.getDomainSegment(p);

        hrleVectorType<hrleIndexType, D> startVector =
            (p == 0) ? grid.getMinGridPoint()
                : newDomain.getSegmentation()[p - 1];

        hrleVectorType<hrleIndexType, D> endVector =
            (p != static_cast<int>(newDomain.getNumberOfSegments() - 1))
                ? newDomain.getSegmentation()[p]
                : grid.incrementIndices(grid.getMaxGridPoint()); 


        //Fast Marching
        for (hrleSparseStarIterator<typename lsDomain<T, D>::DomainType>
          neighborIt(domain, startVector);
          neighborIt.getIndices() < endVector; neighborIt.next()) {

          auto &centerIt = neighborIt.getCenter();

          //mark all active grid points as accepted
          if((activePoints.find(centerIt.getStartIndices()) != activePoints.end())){
            //add them to the new lvl set
            domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                centerIt.getValue());  

          }else{    

            //calculate distance at the current grid point to the interface using the eikonal equation 
            T dist = calcDist(neighborIt, (centerIt.getValue() < 0.));
        
            //check if a distance could be calculated for the current point
            //only points directly on the interface can have a distance of 0 which is not possible here
            if(dist != 0.){

              domainSegment.insertNextDefinedPoint(neighborIt.getIndices(), dist);

            }else{

              domainSegment.insertNextUndefinedPoint(neighborIt.getIndices(), 
              (centerIt.getValue()<0.) ? lsDomain<T, D>::NEG_VALUE : lsDomain<T, D>::POS_VALUE);

            }
            
          }
        }
      } 
      newDomain.finalize();
      levelSet->deepCopy(newlsDomain);
    }


    //TODO: constant is stupid change to dynamic value e.g stop condition or fixed amount of steps
    T width = 3.;
    levelSet->getDomain().segment();
    levelSet->finalize(width);

  }


  private:

  T calcDist(hrleSparseStarIterator<typename lsDomain<T, D>::DomainType>& starStencil, bool inside){

    //TODO: write eikonal equation somwhere

    T stencilMin[D];

    int numUndefined = 0;


    //find the maximum values in the stencil
    for(int i = 0; i < D; i++){

      T pos = starStencil.getNeighbor(i).getValue();

      T neg = starStencil.getNeighbor(i+D).getValue();

      //change the sign to make fast marching on the inside correct
      if(inside){
        pos = pos * -1.;
        neg = neg * -1.;
      }

      //check if the current point in the stencil is defined (not +/- inf)
      if( pos !=  lsDomain<T, D>::POS_VALUE && 
          pos !=  lsDomain<T, D>::NEG_VALUE){

        if(neg !=  lsDomain<T, D>::POS_VALUE && 
           neg !=  lsDomain<T, D>::NEG_VALUE){
            if(pos < neg){
              stencilMin[i] = pos;
            }else{
              stencilMin[i] = neg;
            }
        }else{
          stencilMin[i] = pos;
        }

      }else if(neg !=  lsDomain<T, D>::POS_VALUE && 
               neg !=  lsDomain<T, D>::NEG_VALUE){
        stencilMin[i] = neg;
      }else{
        numUndefined++;
        stencilMin[i] = 0.;
      }

    }

    //all points in the stencil are undefined return 0. (infinity)
    if(numUndefined == D)
      return 0.;

    T sol = 0.;
    
    for(int run = D; run > 0; run--){
      
      T b = 0.;
      T c = (T)run*gridDeltaSqared;

      //calculate upwind difference
      for(int i = 0; i < D; i++){
        b += stencilMin[i];
        c -= (T)run * (stencilMin[i] * stencilMin[i]); 
      }

      T discriminant = (b*b) + c;

      //discriminant is negative must perform a lower dimensional update
      if(discriminant < 0.){
        //remove max value
        T max = 0.;
        int maxIndex = 0;
        for(int i = 0; i < D; i++){
          if(stencilMin[i] > max){
            max = stencilMin[i];
            maxIndex = i;
          }
        }
        stencilMin[maxIndex] = 0.;

        //Test
        if(run == 2){
          T val = 0.;
          for(int i = 0; i < D; i++){
            val += stencilMin[i];
          }
          return (val + 1);
        }
      }else{

        //the distance has to increase and most of the time be greater 0 so we only need the bigger solution
        sol = (b + std::sqrt(discriminant))/(T)run;
        break;
      }
    }    

    if(inside){
      return -sol;
    }else{
      return sol;
    }
    
  }

};



#endif // LS_EIKONAL_EXPAND_HPP