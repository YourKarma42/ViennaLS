#include <iostream>

#include <lsDomain.hpp>

#include <hrleSparseBoxIterator.hpp>

#include <chrono>



template <class T, int D> class silvacoLikeFlagging{

    private:

    typedef typename lsDomain<T, D>::DomainType hrleDomainType;

    lsDomain<T, D> *levelSet = nullptr;

    T cosAngleTreshold = 0.;

    std::vector<hrleVectorType<hrleIndexType, D>> combinations;


    std::unordered_map<hrleVectorType<hrleIndexType,D>, T, typename hrleVectorType<hrleIndexType, D>::hash> flaggedCells;

    public:

    silvacoLikeFlagging(lsDomain<T, D> &passedLevelSet, T passedAngle)
      : levelSet(&passedLevelSet){
      
      cosAngleTreshold = std::cos(passedAngle);
      //TODO: remove rethink
      for(int i = 0; i < D; i ++){
      
        hrleVectorType<hrleIndexType, D> posUnit(0);
        hrleVectorType<hrleIndexType, D> negUnit(0);

        //hrleVectorType<hrleIndexType, D> test = neighborIterator.getIndices();

        int first_pos = i;
        int second_pos = (i+1)%D;
        

        posUnit[first_pos] = 1;
        negUnit[first_pos] = -1;

        combinations.push_back(posUnit);
        combinations.push_back(negUnit);


        posUnit[second_pos] = 1;
        negUnit[second_pos] = 1;

        combinations.push_back(posUnit);
        combinations.push_back(negUnit);

        posUnit[second_pos] = -1;
        negUnit[second_pos] = -1;

        combinations.push_back(posUnit);
        combinations.push_back(negUnit);

      }

      combinations.push_back(hrleVectorType<hrleIndexType, D>(1,1,1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(-1,-1,-1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(-1,-1,1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(1,-1,-1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(-1,1,-1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(-1,1,1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(1,-1,1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(1,1,-1));

      flaggedCells.reserve(levelSet->getNumberOfPoints());

    }

    void createFlagsOutput(std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & activePoints){

        std::vector<T> flags;

        flags.reserve(levelSet->getNumberOfPoints());

        for (hrleConstSparseIterator<hrleDomainType> it(levelSet->getDomain());
          !it.isFinished(); ++it) {

          if (!it.isDefined() || (activePoints.find(it.getStartIndices()) == activePoints.end())) {
              continue;
          }

          if (flaggedCells.find(it.getStartIndices()) == flaggedCells.end()) {
            flags.push_back(0);
          }else{
            flags.push_back(1);
          }
        
        }

        levelSet->getPointData().insertNextScalarData(flags, "Flags Silvaco");
    }


    void apply(std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & activePoints){

      //TODO: remove
      std::vector<std::array<T,D>> outputNormals;
      
      //hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType> neighborIt(levelSet->getDomain());

      std::unordered_map<hrleVectorType<hrleIndexType,D>, std::array<T, D>, typename hrleVectorType<hrleIndexType, D>::hash> normals;

      normals.reserve(activePoints.size());


      auto grid = levelSet->getGrid();

      typename lsDomain<T, D>::DomainType &domain = levelSet->getDomain();

      std::vector<std::unordered_map<hrleVectorType<hrleIndexType, D>, T, typename hrleVectorType<hrleIndexType, D>::hash>> flagsReserve(
      levelSet->getNumberOfSegments());


      double pointsPerSegment =
      double(2 * levelSet->getDomain().getNumberOfPoints()) /
      double(levelSet->getLevelSetWidth());


#pragma omp parallel num_threads((levelSet)->getNumberOfSegments())
        {
            int p = 0;
#ifdef _OPENMP
            p = omp_get_thread_num();
#endif

        hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType> neighborIt(domain);

        std::unordered_map<hrleVectorType<hrleIndexType, D>, std::array<T, D>, typename hrleVectorType<hrleIndexType, D>::hash> normalsSegment;
        normalsSegment.reserve(pointsPerSegment);

        std::unordered_map<hrleVectorType<hrleIndexType, D>, T, typename hrleVectorType<hrleIndexType, D>::hash> &flagsSegment = 
        flagsReserve[p];
        flagsSegment.reserve(pointsPerSegment);

        hrleVectorType<hrleIndexType, D> startVector =
        (p == 0) ? grid.getMinGridPoint()
            : domain.getSegmentation()[p - 1];

        hrleVectorType<hrleIndexType, D> endVector =
        (p != static_cast<int>(domain.getNumberOfSegments() - 1))
            ? domain.getSegmentation()[p]
            : grid.incrementIndices(grid.getMaxGridPoint());



        for(hrleSparseIterator<typename lsDomain<T, D>::DomainType> it(
            domain, startVector);
            it.getStartIndices() < endVector; ++it){

          if (!it.isDefined() || (activePoints.find(it.getStartIndices()) == activePoints.end())) {
            continue;
          }

          neighborIt.goToIndices(it.getStartIndices());

          std::array<T, D> n;

          for (int i = 0; i < D; i++) {
          
            T pos1 = neighborIt.getNeighbor(i).getValue();

            T neg1 = neighborIt.getNeighbor(i + D).getValue();

            n[i] = (pos1 - neg1) * 0.5;
          
          }

          //push normals into a hash map
          normalsSegment[neighborIt.getCenter().getStartIndices()] = n;

          //outputNormals.push_back(n);

          //add to list of normal vectors;

        }
#pragma omp critical
{
        normals.insert(normalsSegment.begin(), normalsSegment.end());
}

#pragma omp barrier

//#pragma omp single
//{


//}

        for(hrleSparseIterator<typename lsDomain<T, D>::DomainType> it(
            domain, startVector);
            it.getStartIndices() < endVector; ++it){

          if (!it.isDefined() || (activePoints.find(it.getStartIndices()) == activePoints.end())) {
              continue;
          } 

          //neighborIterator.goToIndicesSequential(it.getStartIndices());

          std::array<T, D> centerNormal = normals[it.getStartIndices()];

          T normCenterNormal = 0.;

          for(int j= 0; j < D; j++)
            normCenterNormal += centerNormal[j] * centerNormal[j];


          normCenterNormal = std::sqrt(normCenterNormal);

          for(auto dir: combinations){

            if(normals.find(it.getStartIndices() + dir) != normals.end()){
            
              std::array<T, D> currentNormal = normals[it.getStartIndices() + dir];

              T skp = 0.;
              T normCurrentNormal = 0.;

              //calculate scalar product
              for(int j = 0; j < D; j++){

                skp += currentNormal[j]*centerNormal[j];

                normCurrentNormal += currentNormal[j]*currentNormal[j]; 

              }

              normCurrentNormal = std::sqrt(normCurrentNormal);
                            
              if((normCurrentNormal * normCenterNormal * cosAngleTreshold - skp) >= 0.){
                
                flagsSegment[it.getStartIndices()] = 1;

                for(auto dir1: combinations){
                  flagsSegment[it.getStartIndices() + dir1] = 1;
                }

                break;
              }
            }
          }
        }
      }

      for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i) 
        flaggedCells.insert(flagsReserve[i].begin(), flagsReserve[i].end());
      //levelSet->getPointData().insertNextVectorData(outputNormals, "normals");

    }

};