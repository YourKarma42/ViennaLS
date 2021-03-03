#include <iostream>

#include <lsDomain.hpp>

#include <hrleSparseBoxIterator.hpp>

#include <chrono>



template <class T, int D> class silvacoLikeFlagging{

    private:

    typedef typename lsDomain<T, D>::DomainType hrleDomainType;

    lsSmartPointer<lsDomain<double, D>> levelSet = nullptr;

    T cosAngleTreshold = 0.;

    std::vector<hrleVectorType<hrleIndexType, D>> combinations;


    std::vector<T> flaggedCells;

    std::vector<std::array<T, D>> normals;

    public:

    silvacoLikeFlagging(lsSmartPointer<lsDomain<double, D>> passedLevelSet, T passedAngle)
      : levelSet(passedLevelSet){
      
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

    }

    void createFlagsOutput(){

        // insert into pointData of levelSet
        auto &pointData = levelSet->getPointData();
        auto vectorDataPointer = pointData.getScalarData("Angle");
        // if it does not exist, insert new normals vector
        if (vectorDataPointer == nullptr) {
            pointData.insertNextScalarData(flaggedCells, "Angle");
        } else {
        // if it does exist, just swap the old with the new values
            *vectorDataPointer = std::move(flaggedCells);
        }


        std::vector<std::array<T,D>> outputNormals;

        for (hrleConstSparseIterator<hrleDomainType> it(levelSet->getDomain());
          !it.isFinished(); ++it) {
          if (!it.isDefined() || std::abs(it.getValue()) > 0.5) {
            continue;
          }
          outputNormals.push_back(normals[it.getPointId()]);
        }

        pointData.insertNextVectorData(outputNormals, "normals");


    }


    void apply(){

      //clear results from previous run
      normals.clear();
      flaggedCells.clear();

      //std::unordered_map<hrleVectorType<hrleIndexType,D>, std::array<T, D>, typename hrleVectorType<hrleIndexType, D>::hash> normals;

      std::vector<std::vector<std::array<T, D>>> normalsVector(
      levelSet->getNumberOfSegments());


      normals.reserve(levelSet->getNumberOfPoints());
      flaggedCells.reserve(levelSet->getNumberOfPoints());


      auto grid = levelSet->getGrid();

      typename lsDomain<T, D>::DomainType &domain = levelSet->getDomain();

      std::vector<std::vector<T>> flagsReserve(
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
        
        std::array<T,D> zeroVector;

        for(int i = 0; i<D; i++){
          zeroVector[i] = 0.;
        }

        hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType> neighborIt(domain);

        auto &normalsSegment = normalsVector[p];
        normalsSegment.reserve(pointsPerSegment);

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

          if (!it.isDefined()){
              continue;
          }else if(std::abs(it.getValue()) >= 0.5) {
              normalsSegment.push_back(zeroVector);
              continue;
          }

          neighborIt.goToIndices(it.getStartIndices());

          std::array<T, D> n;

          T norm = 0.;

          for (int i = 0; i < D; i++) {
          
            T pos1 = neighborIt.getNeighbor(i).getValue();

            T neg1 = neighborIt.getNeighbor(i + D).getValue();

            n[i] = (pos1 - neg1) * 0.5;

            norm += n[i]*n[i];
          
          }

          norm = std::sqrt(norm);

          for(int j = 0; j < D; j++)
            n[j] = n[j]/norm;
                 

          //push normals into a hash map
          normalsSegment.push_back(n);

        }
      }


      for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i) 
        normals.insert(normals.end(), normalsVector[i].begin(), normalsVector[i].end());


#pragma omp parallel num_threads((levelSet)->getNumberOfSegments())
        {
            int p = 0;
#ifdef _OPENMP
            p = omp_get_thread_num();
#endif

        std::array<T,D> zeroVector;

        for(int i = 0; i<D; i++){
          zeroVector[i] = 0.;
        }

        std::vector<T> &flagsSegment = flagsReserve[p];
        flagsSegment.reserve(pointsPerSegment);

        hrleVectorType<hrleIndexType, D> startVector =
        (p == 0) ? grid.getMinGridPoint()
            : domain.getSegmentation()[p - 1];

        hrleVectorType<hrleIndexType, D> endVector =
        (p != static_cast<int>(domain.getNumberOfSegments() - 1))
            ? domain.getSegmentation()[p]
            : grid.incrementIndices(grid.getMaxGridPoint());

        hrleSparseBoxIterator<hrleDomain<T, D>> boxIterator(domain, 1);

        //hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType> boxIterator(domain);


        for(hrleSparseIterator<typename lsDomain<T, D>::DomainType> it(
            domain, startVector);
            it.getStartIndices() < endVector; ++it){

          if (!it.isDefined() || std::abs(it.getValue()) >= 0.5) {
              continue;
          } 

          std::array<T, D> centerNormal = normals[it.getPointId()];

          boxIterator.goToIndices(it.getStartIndices());
          //std::cout << it.getValue() << std::endl;
          bool flag = false;

         for(unsigned dir = 0; dir < 27; dir++){
          // for(auto dir: combinations){
            //std::cout << boxIterator.getNeighbor(dir).getValue() << std::endl;
            std::array<T, D> currentNormal;
            if(boxIterator.getNeighbor(dir).isDefined()){
              std::array<T, D> currentNormal = normals[boxIterator.getNeighbor(dir).getPointId()];
            }else{
              std::cout << "bad" << std::endl;
              continue;
            }

            if(currentNormal != zeroVector){
            
              T skp = 0.;

              //calculate scalar product
              for(int j = 0; j < D; j++){
                skp += currentNormal[j]*centerNormal[j];
              }

              //vectors are normlized so skp = cos(alpha)          
              if((skp - cosAngleTreshold) >= 0.){
                  flag = true;
                  break;
              }

            }
          }

          if(flag){
            flagsSegment.push_back(1);
          }else{
            flagsSegment.push_back(0);
          }
        }
      }
     
      for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i) 
          flaggedCells.insert(flaggedCells.end(),flagsReserve[i].begin(), flagsReserve[i].end());

    }

};