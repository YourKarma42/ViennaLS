#ifndef LS_EIKONAL_EXPAND_TEST_HPP
#define LS_EIKONAL_EXPAND_TEST_HPP


#include <lsPreCompileMacros.hpp>

#include <hrleSparseIterator.hpp>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>

#include <queue>

#include <lsToMesh.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

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

  //These values are set in the applay function
  T gridDeltaSqared = 0.;
  T gridDelta = 0.;

  int width = 0.;
  int order = 0.;


public:
  lsEikonalExpandTest() {}
  //the band should have at least width 3 for finite differences to work properly
  lsEikonalExpandTest(lsSmartPointer<lsDomain<T, D>> passedLevelSet, 
                      int passedWidth = 3)
      : levelSet(passedLevelSet), width(passedWidth) {
      width = width*2;
      order = width + 2;
  }

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  //TODO: check why lambda dostn work
  //auto myCompare = [](auto &a, auto &b){return (a.first > b.first);}
  struct myCompare{
    bool operator()(std::pair<T, hrleVectorType<hrleIndexType, D>> &a, std::pair<T, hrleVectorType<hrleIndexType, D>> &b){
        return (a.first > b.first);
    }
  };
 
  //velocity extention
  void apply(std::vector<double> &velocities) {

    //the velocity is dependend on the grid resolution
    gridDeltaSqared = levelSet->getGrid().getGridDelta() * levelSet->getGrid().getGridDelta();
    gridDelta = levelSet->getGrid().getGridDelta();

    std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>> container;
    container.reserve(levelSet->getDomain().getNumberOfPoints());

    std::priority_queue <std::pair<T, hrleVectorType<hrleIndexType, D>>, 
                         std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>>, myCompare> 
                          minHeap(myCompare(), std::move(container));

    //collect all grid points needed for velocity extention
      for (hrleSparseStarIterator<hrleDomainType> it(levelSet->getDomain());
        !it.isFinished(); ++it) {

        if(!it.getCenter().isDefined() || std::abs(velocities[it.getCenter().getPointId()]) != 1000)
          continue;

        std::array<T, D> stencilMin;

        //find the min values in the stencil
        for(int i = 0; i < D; i++){

          T pos = it.getNeighbor(i).getValue(); // * gridDelta;

          T neg = it.getNeighbor(i+D).getValue(); // * gridDelta;

          //change the sign to make fast marching on the inside correct
          if(it.getCenter().getValue() < 0.){
            pos = pos * -1.;
            neg = neg * -1.;
          }

          if(pos < neg){
            stencilMin[i] = pos;
          }else{
            stencilMin[i] = neg;
          }

        }

        T dist = solveEikonal(stencilMin, D);
      
        minHeap.push(std::make_pair(std::abs(dist), it.getIndices()));
   
     }

     hrleSparseStarIterator<typename lsDomain<T, D>::DomainType> neighborStarIterator(levelSet->getDomain());

     while(!minHeap.empty()){

        auto currentPoint = minHeap.top();

        minHeap.pop();

        neighborStarIterator.goToIndices(currentPoint.second);

        std::array<T, D> stencilMin;

        //find the min values in the stencil
        for(int i = 0; i < D; i++){

          T pos = neighborStarIterator.getNeighbor(i).getValue(); // * gridDelta;

          T neg = neighborStarIterator.getNeighbor(i+D).getValue(); // * gridDelta;

          //change the sign to make fast marching on the inside correct
          if(neighborStarIterator.getCenter().getValue() < 0.){
            pos = pos * -1.;
            neg = neg * -1.;
          }

          if(pos < neg){
            stencilMin[i] = velocities[neighborStarIterator.getNeighbor(i).getPointId()];
          }else{
            stencilMin[i] = velocities[neighborStarIterator.getNeighbor(i+D).getPointId()];
          }

        }

        T dist = solveEikonal(stencilMin, D);

        velocities[neighborStarIterator.getCenter().getPointId()] = dist;

        /*if(neighborStarIterator.getCenter().getValue() < 0.){
          velocities[neighborStarIterator.getCenter().getPointId()] = -dist;
        }else{
          velocities[neighborStarIterator.getCenter().getPointId()] = dist;
        }*/

                  
     }
  }


  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsEikonalExpand.")
          .print();
      return;
    }

    //LS values are normalized between 0 and 1
    gridDelta = 1;
    gridDeltaSqared = 1;

    //currently Eikonal expand is only single core so the ls has to be contained in one segment
    levelSet->getDomain().desegment();

    accepted.clear();

    //used to keep track of the Interface
    T largeValue = 1000.;

    //minimum heap for FMM
    std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>> container;
    container.reserve(levelSet->getDomain().getNumberOfPoints() * 10);
    //std::priority_queue <std::pair<T, hrleVectorType<hrleIndexType, D>>, std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>>, myCompare> minHeap;

    std::priority_queue <std::pair<T, hrleVectorType<hrleIndexType, D>>, 
                         std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>>, 
                         myCompare> 
                            minHeap(myCompare(), std::move(container));

    //std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>> initialValues;

    //expand the domain for FMM
    for(int runs = 0; runs < order; runs++){

        auto &grid = levelSet->getGrid();
        auto newlsDomain = lsSmartPointer<lsDomain<T, D>>::New(grid, lsNormalizations::EUCLID);
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
                    //interface.insert(it.getIndices());
                }

            }else{

            int undefined = 0;

            for(int i = 0; i < D; i++){

              T pos = std::abs(it.getNeighbor(i).getValue());

              T neg = std::abs(it.getNeighbor(i+D).getValue());

              if(pos !=  lsDomain<T, D>::POS_VALUE || neg !=  lsDomain<T, D>::POS_VALUE)
                break;
              
              undefined++;
            }

            if(undefined < D){
                newDomain.getDomainSegment(0).insertNextDefinedPoint(it.getIndices(),
                                            (it.getCenter().getValue()<0.) ? -largeValue : largeValue);                               
            }else{
                newDomain.getDomainSegment(0).insertNextUndefinedPoint(it.getIndices(), 
                    (it.getCenter().getValue()<0.) ? lsDomain<T, D>::NEG_VALUE : lsDomain<T, D>::POS_VALUE);
            }
          }
        }
        if(runs == 0)
          largeValue = largeValue+100;
        levelSet->deepCopy(newlsDomain);     
    }

    auto narrowband0 = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<T, D>(levelSet, narrowband0, true, true, 2*largeValue+10).apply();

    lsVTKWriter(narrowband0, lsFileFormatEnum::VTU , "afterExpansion" ).apply();

    std::vector<double> tmp;
    //initialize the Domain for FMM
    for (hrleSparseIterator<hrleDomainType> it(levelSet->getDomain());
        !it.isFinished(); ++it) {

        if(!it.isDefined())
          continue;

        if(std::abs(it.getValue()) <= largeValue-100){
            accepted.push_back(2);
            tmp.push_back(2);

        }else{
            accepted.push_back(0); 
            tmp.push_back(0); 
        }                           
    }


    //initialize the narrowband for FMM
    for (hrleSparseStarIterator<hrleDomainType> it(levelSet->getDomain());
        !it.isFinished(); ++it) {

        if(!it.getCenter().isDefined() || std::abs(it.getCenter().getValue()) != (largeValue-100))
          continue;
       

        //TODO: calc dist only accepted
        T dist = calcDistAccepted(it, (it.getCenter().getValue() < 0.));   
       
        if (dist != 0){

            T &currentVal = it.getCenter().getValue();
          	currentVal = dist;

            minHeap.push(std::make_pair(std::abs(dist), it.getIndices()));
            accepted[it.getCenter().getPointId()] = 1;
            tmp[it.getCenter().getPointId()] = 1;
        }

    }

/*
    auto narrowband1 = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<T, D>(levelSet, narrowband1, true, true, largeValue).apply();

    narrowband1->insertNextScalarData(tmp, "accepted");

    lsVTKWriter(narrowband1, lsFileFormatEnum::VTU , "beforMarch" ).apply();
*/

    march(minHeap);

/*
    auto narrowband = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<T, D>(levelSet, narrowband, true, true, 1000).apply();

    narrowband->insertNextScalarData(tmp, "accepted");

    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "afterMarch" ).apply();
*/

    levelSet->getDomain().segment();

    levelSet->finalize(width);
  }

  private:

  //0 ... far
  //1 ... considered
  //2 ... accepted
  std::vector<int> accepted;


  void march(std::priority_queue <std::pair<T, hrleVectorType<hrleIndexType, D>>, 
                         std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>>, 
                         myCompare> & minHeap ){

      hrleSparseStarIterator<typename lsDomain<T, D>::DomainType> neighborStarIterator(levelSet->getDomain());

      hrleSparseStarIterator<typename lsDomain<T, D>::DomainType> starItEikonal(levelSet->getDomain());

      while(!minHeap.empty()){

          auto currentPoint = minHeap.top();

          minHeap.pop();

          neighborStarIterator.goToIndices(currentPoint.second);

          //is the current value already accepted
          if(accepted[neighborStarIterator.getCenter().getPointId()] <= 1){
            

              //break the loop if the narrowband has expanded far enough for application
              if(currentPoint.first > width)
                break;

              accepted[neighborStarIterator.getCenter().getPointId()] = 2;

              //update neighbours
              for(int i = 0; i < 2*D; i++){

                  auto currentNeighbor = neighborStarIterator.getNeighbor(i);

                  //if curretn neighbour is undefined skip (points at the border of the domain)
                  if(currentNeighbor.isDefined()){

                      //if current point is accepted skip it
                      if(accepted[currentNeighbor.getPointId()] <= 1){

                          T &currentValue = currentNeighbor.getValue();

                          starItEikonal.goToIndices(currentNeighbor.getOffsetIndices());
                          
                          //send iterator
                          T dist = calcDistFMM(starItEikonal, (currentNeighbor.getValue() < 0.));   


                          if(std::abs(dist) < std::abs(currentValue)){

                              currentValue = dist;
                              minHeap.push(std::make_pair(std::abs(dist), currentNeighbor.getOffsetIndices()));                           
                          }                 
                          //TODO ERROR?
                          accepted[currentNeighbor.getPointId()] = 1;
                          
                      }

                  }
              }

          }
      }


  }

  T calcDistAccepted(hrleSparseStarIterator<typename lsDomain<T, D>::DomainType>& starStencil, bool inside){

    std::array<T,D> stencilMin;

    int dim = 0;
    int undefined = 0;

    //find the maximum values in the stencil
    for(int i = 0; i < D; i++){

      T pos = starStencil.getNeighbor(i).getValue(); // * gridDelta;

      T neg = starStencil.getNeighbor(i+D).getValue(); // * gridDelta;


      //change the sign to make fast marching on the inside correct
      if(inside){
        pos = pos * -1.;
        neg = neg * -1.;
      }

      //Take care of undefined points (not in accepted array)

      int posAccepted = 0;
      if(starStencil.getNeighbor(i).isDefined() && pos != lsDomain<T, D>::POS_VALUE){
          posAccepted = accepted[starStencil.getNeighbor(i).getPointId()];
      }


      int negAccepted = 0;
      if(starStencil.getNeighbor(i+D).isDefined() && neg != lsDomain<T, D>::POS_VALUE){
          negAccepted = accepted[starStencil.getNeighbor(i+D).getPointId()];
      }

      if(posAccepted == 2){        
        if(negAccepted == 2){
            if(pos < neg){
              stencilMin[i] = pos;
            }else{
              stencilMin[i] = neg;
            }
        }else{
          stencilMin[i] = pos;
        }
        dim++;
      }else if(negAccepted == 2){
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

    //solve Eikonal equation
    T sol = solveEikonal(stencilMin, dim);

    if(inside){
      return -sol; ///gridDelta;
    }else{
      return sol; ///gridDelta;
    }

  }



  T calcDistFMM(hrleSparseStarIterator<typename lsDomain<T, D>::DomainType>& starStencil, bool inside){

    //TODO: write eikonal equation somwhere

    std::array<T,D> stencilMin;

    int dim = 0;
    int undefined = 0;

    //find the maximum values in the stencil
    for(int i = 0; i < D; i++){

      T pos = starStencil.getNeighbor(i).getValue(); // * gridDelta;

      T neg = starStencil.getNeighbor(i+D).getValue(); // * gridDelta;


      //change the sign to make fast marching on the inside correct
      if(inside){
        pos = pos * -1.;
        neg = neg * -1.;
      }

      //Take care of undefined points (not in accepted array)

      int posAccepted = 0;
      if(starStencil.getNeighbor(i).isDefined() && pos != lsDomain<T, D>::POS_VALUE){
          posAccepted = accepted[starStencil.getNeighbor(i).getPointId()];
      }


      int negAccepted = 0;
      if(starStencil.getNeighbor(i+D).isDefined() && neg != lsDomain<T, D>::POS_VALUE){
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

    //solve Eikonal equation
    T sol = solveEikonal(stencilMin, dim);

    if(inside){
      return -sol; ///gridDelta;
    }else{
      return sol; ///gridDelta;
    }
    
  }

  T solveEikonal(std::array<T, D>& stencilMin, int dim){
    T sol = 0.;
    
    for(int run = D; run > 0; run--){
      //perform one dimensional update
      if(dim == 1){
          for(int i = 0; i < D; i++)
              sol += stencilMin[i];
          return sol + gridDelta;
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
        //the distance has to increase so we only need the bigger solution
        return (-b + std::sqrt(q))/(2.*a);
        break;
      }
    }    
  }

};



#endif // LS_EIKONAL_EXPAND_TEST_HPP

