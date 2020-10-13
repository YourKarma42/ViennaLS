#ifndef LS_EXPAND_SPHERE_HPP
#define LS_EXPAND_SPHERE_HPP


#include <lsPreCompileMacros.hpp>

#include <hrleSparseIterator.hpp>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>

#include <unordered_set>

//TODO: expand!
//Important: This implementation uses a BRUTE FORCE approach to calculate the signed distance field.
//To achive better performance a fast sweeping approach should be used.
/// This class expands the level set to a certain numeber of layers by using the Eikonal equation.
/// This expansion is more accurate than using Manhatten expansion, however it has alonger Calculation time.
/// This form of expansion keeps the signed distance porperty of the level set |grad(f)| = 1.

template <class T, int D> class lsExpandSphere {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;

  T gridDelta = 0.;

  T radius = 0.;

  hrleVectorType<T, D> origin;

  //TOD: make call by reference (problem with initialization)
  std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash>  activePoints;


public:
  lsExpandSphere() {}

  lsExpandSphere(lsSmartPointer<lsDomain<T, D>> passedLevelSet, 
      std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash>  passedActivePoints, 
      T passedRadius, hrleVectorType<T, D> passedOrigin)
      : levelSet(passedLevelSet), activePoints(passedActivePoints), radius(passedRadius), origin(passedOrigin){
        //TODO: check if level set is euler normalized and give error if not
      gridDelta = levelSet->getGrid().getGridDelta();
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
            1 + 1.0 / static_cast<double>(runs+1);
        //TODO:probably remove
        //const T limit = (runs + 1) * T(0.5);
        
        auto &grid = levelSet->getGrid();

        auto newlsDomain = lsSmartPointer<lsDomain<T, D>>::New(grid);
        typename lsDomain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
        typename lsDomain<T, D>::DomainType &domain = levelSet->getDomain();

        newDomain.initialize(domain.getNewSegmentation(),
                            domain.getAllocation() * allocationFactor);

        std::cout << "min index " << grid.getMinIndex() << std::endl;

        std::cout << "max index " << grid.getMaxIndex() << std::endl;

//#pragma omp parallel num_threads(newDomain.getNumberOfSegments())
        {
            int p = 0;
//#ifdef _OPENMP
            //p = omp_get_thread_num();
//#endif

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



                int numUndefined = 0;

                auto &centerIt = neighborIt.getCenter();

                for(int i = 0; i < D; i++){

                  T pos = neighborIt.getNeighbor(i).getValue();

                  T neg = neighborIt.getNeighbor(i+D).getValue();

                  //change the sign to make fast marching on the inside correct
                 // if(inside){
                 //   pos = pos * -1.;
                 //   neg = neg * -1.;
                 // }

                  //check if the current point in the stencil is defined (not +/- inf)
                  //TODO:this if is stupid
                  if( pos !=  lsDomain<T, D>::POS_VALUE && 
                      pos !=  lsDomain<T, D>::NEG_VALUE){

                  }else if(neg !=  lsDomain<T, D>::POS_VALUE && 
                          neg !=  lsDomain<T, D>::NEG_VALUE){

                  }else{
                    numUndefined++;
                  }

                }

                //TODO: why was there an error when not using second condition?
                if(numUndefined == D && 
                  (activePoints.find(centerIt.getStartIndices()) == activePoints.end())){

                  domainSegment.insertNextUndefinedPoint(neighborIt.getIndices(), 
                    (centerIt.getValue()<0.) ? lsDomain<T, D>::NEG_VALUE : lsDomain<T, D>::POS_VALUE);
                }else{
                  T pointRadius = 0.;

                  for(int i = 0; i < D; i++){
                      //std::cout << neighborIt.getIndices()[i] << std::endl;
                      pointRadius += ( gridDelta * neighborIt.getIndices()[i] - origin[i]) * ( gridDelta * neighborIt.getIndices()[i] - origin[i]);
                  }

                  pointRadius = std::sqrt(pointRadius);

                  T dist = pointRadius - radius;                       
      
                  domainSegment.insertNextDefinedPoint(neighborIt.getIndices(), dist);

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


};



#endif // LS_EXPAND_SPHERE_HPP