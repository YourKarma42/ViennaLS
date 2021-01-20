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
        //TODO: check if level set is euler normalized and give error if not
      gridDeltaSqared = levelSet->getGrid().getGridDelta() * levelSet->getGrid().getGridDelta();
      gridDelta = levelSet->getGrid().getGridDelta();
      //TODO:discuss if +1 sufficiant
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

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsEikonalExpand.")
          .print();
      return;
    }

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

    std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>> initialValues;

    //expand the domain for FMM
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
        levelSet->deepCopy(newlsDomain);     
    }

    std::vector<double> tmp;
    //initialize the Domain for FMM
    for (hrleSparseIterator<hrleDomainType> it(levelSet->getDomain());
        !it.isFinished(); ++it) {

        if(!it.isDefined())
          continue;

        if(std::abs(it.getValue()) <= gridDelta){
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

        if(!it.getCenter().isDefined() || std::abs(it.getCenter().getValue()) <= gridDelta)
          continue;

        T dist = calcDistFMM(it, (it.getCenter().getValue() < 0.));   


        if ( dist != 0){

            T &currentVal = it.getCenter().getValue();

          	currentVal = dist;

            minHeap.push(std::make_pair(std::abs(dist), it.getIndices()));

        }

    }
/*
    auto narrowband1 = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<T, D>(levelSet, narrowband1, true, true, largeValue + 1).apply();

    narrowband1->insertNextScalarData(tmp, "accepted");

    lsVTKWriter(narrowband1, lsFileFormatEnum::VTU , "beforMarch" ).apply();
*/
    //check min heap
    march(minHeap);
/*
    auto narrowband = lsSmartPointer<lsMesh>::New();
    std::cout << "Extracting narrowband..." << std::endl;
    lsToMesh<T, D>(levelSet, narrowband, true, true, largeValue + 1).apply();

    narrowband->insertNextScalarData(tmp, "accepted");

    lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "afterMarch" ).apply();
*/

    levelSet->getDomain().segment();

    levelSet->finalize(width);
  }


  void applyTest() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsEikonalExpand.")
          .print();
      return;
    }

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

    std::vector<std::pair<T, hrleVectorType<hrleIndexType, D>>> initialValues;

    //initialize the Domain for FMM
    for (hrleSparseIterator<hrleDomainType> it(levelSet->getDomain());
        !it.isFinished(); ++it) {

        if(!it.isDefined())
          continue;

        if(it.getValue() <= gridDelta){
            accepted.push_back(2);

        }else{
            accepted.push_back(0); 
            T &currentVal = it.getValue();
            currentVal = (it.getCenter().getValue()<0.) ? -largeValue : largeValue;
        }                           
    }

    //initialize the narrowband for FMM
    for (hrleSparseStarIterator<hrleDomainType> it(levelSet->getDomain());
        !it.isFinished(); ++it) {

        if(!it.getCenter().isDefined() && it.getCenter().getValue() > gridDelta)
          continue;

        T dist = calcDistFMM(it, (it.getCenter().getValue() < 0.));   

        if ( dist != 0){

            T &currentVal = it.getCenter().getValue();

          	currentVal = dist;

            minHeap.push(std::make_pair(std::abs(dist)), it.getCenter().getStartIndices());

        }
            

    }

    march(minHeap);

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

          //std::cout << neighborStarIterator.getCenter().isDefined() << " " << 
          //neighborStarIterator.getCenter().getValue() << " " << neighborStarIterator.getCenter().getPointId()<< std::endl; 

          //is the current value already accepted
          //if(accepted.at(neighborStarIterator.getCenter().getPointId()) <= 1){
          if(accepted[neighborStarIterator.getCenter().getPointId()] <= 1){
            

              //break the loop if the narrowband has expanded far enough for application
              if(currentPoint.first > width*gridDelta)
                break;

              accepted[neighborStarIterator.getCenter().getPointId()] = 2;

              //update neighbours
              for(int i = 0; i < 2*D; i++){

                  auto currentNeighbor = neighborStarIterator.getNeighbor(i);

                  //if curretn neighbour is undefined skip (points at the border of the domain)
                  if(currentNeighbor.isDefined()){
                  
                      if(accepted[currentNeighbor.getPointId()] <= 1){

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
      }


  }

  T calcDistFMM(hrleSparseStarIterator<typename lsDomain<T, D>::DomainType>& starStencil, bool inside){

    //TODO: write eikonal equation somwhere

    T stencilMin[D];

    int dim = 0;
    int undefined = 0;
/*

      hrleVectorType<hrleIndexType, D> stopVec(-25,16);

    if(starStencil.getCenter().getStartIndices() == stopVec){
        auto narrowband = lsSmartPointer<lsMesh>::New();
        std::cout << "Extracting FMM result..." << std::endl;
        lsToMesh<T, D>(levelSet, narrowband, true, false, 20).apply();

        std::vector<T> acceptedT(accepted.begin(), accepted.end());

        narrowband->insertNextScalarData(acceptedT, "accepted");
        lsVTKWriter(narrowband, lsFileFormatEnum::VTU , "/media/sf_shared/duringFMM" ).apply();

        std::cout << "pause" <<std::endl;
    }

*/
    //find the maximum values in the stencil
    for(int i = 0; i < D; i++){

      T pos = starStencil.getNeighbor(i).getValue();

      T neg = starStencil.getNeighbor(i+D).getValue();


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