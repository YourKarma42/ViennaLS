#ifndef LS_EIKONAL_EXPAND_TEST_HPP
#define LS_EIKONAL_EXPAND_TEST_HPP


#include <lsPreCompileMacros.hpp>

#include <hrleSparseIterator.hpp>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>

#include <queue>

#include <unordered_set>

#include <unordered_map>

//TODO: expand!
//Important: This implementation uses a BRUTE FORCE approach to calculate the signed distance field.
//To achive better performance a fast sweeping approach should be used.
/// This class expands the level set to a certain numeber of layers by using the Eikonal equation.
/// This expansion is more accurate than using Manhatten expansion, however it has alonger Calculation time.
/// This form of expansion keeps the signed distance porperty of the level set |grad(f)| = 1.

template <class T, int D> class lsEikonalExpandTest {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;

  T gridDeltaSqared = 0.;
  T gridDelta = 0.;

  //TOD: make call by reference (problem with initialization)
  std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash>  activePoints;


public:
  lsEikonalExpandTest() {}

  lsEikonalExpandTest(lsSmartPointer<lsDomain<T, D>> passedLevelSet, 
      std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash>  passedActivePoints)
      : levelSet(passedLevelSet), activePoints(passedActivePoints) {
        //TODO: check if level set is euler normalized and give error if not
      gridDeltaSqared = levelSet->getGrid().getGridDelta() * levelSet->getGrid().getGridDelta();
      gridDelta = levelSet->getGrid().getGridDelta();
  }

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  void setActivePoints(std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash>  passedActivePoints) {
    activePoints = passedActivePoints;
  }


  struct myCompare{
    bool operator()(std::pair<T, hrleVectorType<hrleIndexType, D>> &a, std::pair<T, hrleVectorType<hrleIndexType, D>> &b){
        return (a.first > b.first);
    }
  };

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsEikonalExpand.")
          .print();
      return;
    }

    //TODO: strange bug with 0 element when using to big order
    int order = 6;

    //used to keep track of the Interface
    T largeValue = 1000.;

    std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> interface;

    hrleVectorType<hrleIndexType, D> tmpVec(0);
    std::pair<T, hrleVectorType<hrleIndexType, D>> startValueOutside(largeValue, tmpVec);
    std::pair<T, hrleVectorType<hrleIndexType, D>> startValueInside(largeValue, tmpVec);


    //minimum heap for FMM
    std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>> container;
    container.reserve(levelSet->getDomain().getNumberOfPoints() * 10);
    //std::priority_queue <std::pair<T, hrleVectorType<hrleIndexType, D>>, std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>>, myCompare> minHeap;

    std::priority_queue <std::pair<T, hrleVectorType<hrleIndexType, D>>, 
                         std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>>, 
                         myCompare> 
                            minHeap(myCompare(), std::move(container));

    std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>> initialValues;

    //initialize the Domain for FMM
    for (hrleSparseStarIterator<hrleDomainType> it(levelSet->getDomain());
        !it.isFinished(); ++it) {

        if(it.getCenter().isDefined())
          continue;

        T dist = calcDist(it, (it.getCenter().getValue() < 0.));
        if(dist != 0){
            initialValues.push_back(std::make_pair(dist, it.getIndices()));            
        }                            
    }

    for(int runs = 0; runs < order; runs++){

        auto &grid = levelSet->getGrid();
        auto newlsDomain = lsSmartPointer<lsDomain<T, D>>::New(grid);
        typename lsDomain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
        typename lsDomain<T, D>::DomainType &domain = levelSet->getDomain();

        const int allocationFactor =
          1 + 1.0 / static_cast<double>(runs);

        newDomain.initialize(domain.getNewSegmentation(),
                           domain.getAllocation() * allocationFactor);

        //create initialize FMM
        for (hrleSparseStarIterator<hrleDomainType> it(levelSet->getDomain());
            !it.isFinished(); ++it) {

            if(it.getCenter().isDefined()){
                newDomain.getDomainSegment(0).insertNextDefinedPoint(it.getIndices(),
                                            it.getCenter().getValue());
                if(runs == 0){
                    interface.insert(it.getIndices());
                }

            }else{

              //TODO: rethink can problably made faster not always calc dist
              T dist = calcDist(it, (it.getCenter().getValue() < 0.));
                if(dist != 0){                      
                    newDomain.getDomainSegment(0).insertNextDefinedPoint(it.getIndices(),
                                            (it.getCenter().getValue()<0.) ? -largeValue : largeValue);                  
                }else{
                    newDomain.getDomainSegment(0).insertNextUndefinedPoint(it.getIndices(), 
                    (it.getCenter().getValue()<0.) ? lsDomain<T, D>::NEG_VALUE : lsDomain<T, D>::POS_VALUE);
                }
            }
        }
        newDomain.finalize();
        levelSet->deepCopy(newlsDomain);
    }

    //initialize the narrowband for FMM
    for(auto p: initialValues){
        hrleSparseIterator<hrleDomainType> it(levelSet->getDomain());
        it.goToIndices(p.second);
        T &currentVal = it.getValue();

        currentVal = p.first;

        minHeap.push(std::make_pair(std::abs(p.first), p.second));
    }

    T width = order;
    //TODO: carful with segment here multicore!


    levelSet->finalize(width);

    std::cout << "ready to march" << std::endl;


    for (hrleSparseIterator<hrleDomainType> it(levelSet->getDomain());
        !it.isFinished(); ++it) {

        if(it.isDefined()){
            if(interface.find(it.getStartIndices()) != interface.end()){
                accepted.push_back(2);
            }else{
                accepted.push_back(0);     
            }
        }
    }

    std::cout << "test";

    hrleSparseStarIterator<typename lsDomain<T, D>::DomainType> neighborStarIterator(levelSet->getDomain());

    hrleSparseStarIterator<typename lsDomain<T, D>::DomainType> starItEikonal(levelSet->getDomain());

    while(!minHeap.empty()){

        auto currentPoint = minHeap.top();



        minHeap.pop();

        neighborStarIterator.goToIndices(currentPoint.second);

        //is the current value already accepted
        if(accepted[neighborStarIterator.getCenter().getPointId()] <= 1){

            //std::cout << minHeap.size() << std::endl;

            if(currentPoint.first > 3*gridDelta)
              break;

            accepted[neighborStarIterator.getCenter().getPointId()] = 2;

            //update neighbours
            for(int i = 0; i < 2*D; i++){

                auto currentNeighbor = neighborStarIterator.getNeighbor(i);

                int tmpAccepted = 0;

                if(currentNeighbor.isDefined()){
                    tmpAccepted = accepted[currentNeighbor.getPointId()]; 
                }

                //if not accepted
                if(tmpAccepted <= 1){

                    T &currentValue = currentNeighbor.getValue();

                    starItEikonal.goToIndices(currentNeighbor.getOffsetIndices());
                    
                    //send iterator
                    T dist = calcDistFMM(starItEikonal, (currentNeighbor.getValue() < 0.));   


                    if(abs(dist) < abs(currentValue)){

                        currentValue = dist;
                        minHeap.push(std::make_pair(std::abs(dist), currentNeighbor.getOffsetIndices()));                           
                    }                 

                    accepted[currentNeighbor.getPointId()] == 1;
                    

                }
            }

        }
    }

    levelSet->getDomain().segment();

    levelSet->finalize(width);
  }

  private:

  //0 ... far
  //1 ... considered
  //2 ... accepted
  std::vector<int> accepted;


  T calcDistFMM(hrleSparseStarIterator<typename lsDomain<T, D>::DomainType>& starStencil, bool inside){

    //TODO: write eikonal equation somwhere

    T stencilMin[D];

    int dim = 0;
    int undefined = 0;


    //find the maximum values in the stencil
    for(int i = 0; i < D; i++){

      T pos = starStencil.getNeighbor(i).getValue();

      T neg = starStencil.getNeighbor(i+D).getValue();


      //change the sign to make fast marching on the inside correct
      if(inside){
        pos = pos * -1.;
        neg = neg * -1.;
      }


      //STENCIL CAN ACESS MEMORY THAT IS NOT INITIALZED IN accepted

      int posAccepted = 0;
      if(pos != lsDomain<T, D>::POS_VALUE){
          posAccepted = accepted[starStencil.getNeighbor(i).getPointId()];
      }


      int negAccepted = 0;
      if(neg != lsDomain<T, D>::POS_VALUE){
          negAccepted = accepted[starStencil.getNeighbor(i+D).getPointId()];
      }
      



      //check if the current point in the stencil is defined (not +/- inf)


      if( pos !=  lsDomain<T, D>::POS_VALUE && 
          posAccepted == 2){
          //pos !=  lsDomain<T, D>::NEG_VALUE){

        if(neg !=  lsDomain<T, D>::POS_VALUE && 
           //neg !=  lsDomain<T, D>::NEG_VALUE){
            negAccepted == 2){
            if(pos < neg){
              stencilMin[i] = pos;
            }else{
              stencilMin[i] = neg;
            }
        }else{
          stencilMin[i] = pos;
        }
        dim++;
      }else if(neg !=  lsDomain<T, D>::POS_VALUE && 
              negAccepted == 2){
               //neg !=  lsDomain<T, D>::NEG_VALUE){
        stencilMin[i] = neg;
        dim++;
      }else{
        undefined++;
        stencilMin[i] = 0.;
      }

    }

    //all points in the stencil are undefined return 0. (infinity)
    if(undefined == D)
      return 0.;

    T sol = 0.;
    
    for(int run = D; run > 0; run--){
      //perform one dimensional update
      if(dim == 1){
          for(int i = 0; i < D; i++)
              sol += stencilMin[i];
          sol += gridDelta;
          break;
      }

      T sumT = 0.;
      T sumTsq = 0.;

      for(int i = 0; i < D; i++){
        sumT += stencilMin[i];
        sumTsq += stencilMin[i]*stencilMin[i];
      }

      T a = T(dim);
      T b = -2*sumT;
      T c = sumTsq - gridDeltaSqared;
      T q = b*b - 4*a*c;
      
      //discriminant is negative must perform a lower dimensional update
      if(q < 0.){
        //remove max value
        T max = 0.;
        int maxIndex = 0;
        for(int i = 0; i < D; i++){
          if(stencilMin[i] > max){
            max = stencilMin[i];
            maxIndex = i;
          }
        }
        dim--;
        stencilMin[maxIndex] = 0.;

      }else{

        //the distance has to increase and most of the time be greater 0 so we only need the bigger solution
        sol = (-b + std::sqrt(q))/(2.*a);
        break;
      }
    }    

    if(inside){
      return -sol;
    }else{
      return sol;
    }
    
  }

    T calcDist(hrleSparseStarIterator<typename lsDomain<T, D>::DomainType>& starStencil, bool inside){

    //TODO: write eikonal equation somwhere

    T stencilMin[D];

    int dim = 0;
    int undefined = 0;


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
        dim++;
      }else if(neg !=  lsDomain<T, D>::POS_VALUE && 
               neg !=  lsDomain<T, D>::NEG_VALUE){
        stencilMin[i] = neg;
        dim++;
      }else{
        undefined++;
        stencilMin[i] = 0.;
      }

    }

    //all points in the stencil are undefined return 0. (infinity)
    if(undefined == D)
      return 0.;

    T sol = 0.;
    
    for(int run = D; run > 0; run--){
      //perform one dimensional update
      if(dim == 1){
          for(int i = 0; i < D; i++)
              sol += stencilMin[i];
          sol += gridDelta;
          break;
      }

      T sumT = 0.;
      T sumTsq = 0.;

      for(int i = 0; i < D; i++){
        sumT += stencilMin[i];
        sumTsq += stencilMin[i]*stencilMin[i];
      }

      T a = T(dim);
      T b = -2*sumT;
      T c = sumTsq - gridDeltaSqared;
      T q = b*b - 4*a*c;
      
      //discriminant is negative must perform a lower dimensional update
      if(q < 0.){
        //remove max value
        T max = 0.;
        int maxIndex = 0;
        for(int i = 0; i < D; i++){
          if(stencilMin[i] > max){
            max = stencilMin[i];
            maxIndex = i;
          }
        }
        dim--;
        stencilMin[maxIndex] = 0.;

      }else{

        //the distance has to increase and most of the time be greater 0 so we only need the bigger solution
        sol = (-b + std::sqrt(q))/(2.*a);
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



#endif // LS_EIKONAL_EXPAND_TEST_HPP