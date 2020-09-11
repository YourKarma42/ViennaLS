#include <iostream>

#include <lsDomain.hpp>

#include <hrleSparseBoxIterator.hpp>

#include <chrono>

#include <../dev/derivatives.hpp>



template <class T, int D> class myFlagging{

    private:

    typedef typename lsDomain<T, D>::DomainType hrleDomainType;

    lsDomain<T, D> *levelSet = nullptr;

    T flatBoundary = 0.;

    std::unordered_map<hrleVectorType<hrleIndexType,D>, T, typename hrleVectorType<hrleIndexType, D>::hash> flaggedCells;

    std::vector<T> curveOutput;

    public:

    myFlagging(lsDomain<T, D> &passedLevelSet, T passedBoundary)
      : levelSet(&passedLevelSet), flatBoundary(passedBoundary){
        //flags.reserve(levelSet->getNumberOfPoints());

        curveOutput.reserve(levelSet->getNumberOfPoints());

    }

    void createFlagsOutput(std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & activePoints){

        std::vector<T> flags;

        flaggedCells.reserve(levelSet->getNumberOfPoints());   

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

        levelSet->getPointData().insertNextScalarData(flags, "Flags");

        //levelSet->getPointData().insertNextScalarData(curveOutput, "curvature");
    }


    void apply(std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> & activePoints){
        

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


            curvaturGeneralFormula<T, D> curvatureCalculator(levelSet->getGrid().getGridDelta());

            hrleSparseBoxIterator<hrleDomain<T, D>> boxIterator(levelSet->getDomain(), 1);

            //curvaturGeneralFormulaBigStencilBias <T, D> curvatureCalculator(1.0);

            //hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType> neighborIt(domain);

            //curvaturShapeDerivatives1<T, D> curvatureCalculator(levelSet->getGrid().getGridDelta());


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

            //for(hrleConstSparseIterator<hrleDomainType> it(levelSet->getDomain());
             //   !it.isFinished(); ++it){

                if (!it.isDefined() || (activePoints.find(it.getStartIndices()) == activePoints.end())) {
                    continue;
                } 

                //neighborIt.goToIndices(it.getStartIndices());

                //T curve = curvatureCalculator(neighborIt);

                boxIterator.goToIndicesSequential(it.getStartIndices());

                T curve = curvatureCalculator(boxIterator);

                //curveOutput.push_back(curve);

                if(std::abs(curve) > flatBoundary){
                    //flagsSegment.insert({it.getStartIndices(), 1});
                    flagsSegment[it.getStartIndices()] = 1;
                }
                
            } 

        }

        for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i) 
            flaggedCells.insert(flagsReserve[i].begin(), flagsReserve[i].end());
      
   // }

    }

};