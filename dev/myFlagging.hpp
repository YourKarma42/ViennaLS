#include <iostream>

#include <lsDomain.hpp>

#include <hrleSparseBoxIterator.hpp>

#include <chrono>

#include <../dev/derivatives.hpp>
#include <../dev/hrleTestIterator.hpp>



template <class T, int D> class myFlagging{

    private:

    typedef typename lsDomain<T, D>::DomainType hrleDomainType;

    lsSmartPointer<lsDomain<double, D>> levelSet = nullptr;

    T flatBoundary = 0.;

    std::vector<T> flaggedCells;

    std::vector<T> curveOutput;

    public:

    myFlagging(lsSmartPointer<lsDomain<double, D>> passedLevelSet, T passedBoundary)
      : levelSet(passedLevelSet), flatBoundary(passedBoundary){
        //flags.reserve(levelSet->getNumberOfPoints());

        //curveOutput.reserve(levelSet->getNumberOfPoints());

    }

    void createFlagsOutput(int method){

        std::string name; 

        if(method == 0){
            name = "Shape";
            std::cout << "Shape" << std::endl;
        }else{
            name = "General";
            std::cout << "General" << std::endl;
        }

        // insert into pointData of levelSet
        auto &pointData = levelSet->getPointData();
        auto vectorDataPointer = pointData.getScalarData(name);
        // if it does not exist, insert new normals vector
        if (vectorDataPointer == nullptr) {
            pointData.insertNextScalarData(flaggedCells, name);
        } else {
        // if it does exist, just swap the old with the new values
            *vectorDataPointer = std::move(flaggedCells);
        }

    }

    //0 shape operator; 1 general formula
    double apply(int method){

        
        auto start = std::chrono::high_resolution_clock::now(); 
        auto stop = std::chrono::high_resolution_clock::now(); 

        flaggedCells.clear();
        flaggedCells.reserve(levelSet->getNumberOfPoints());
        

        auto grid = levelSet->getGrid();

        typename lsDomain<T, D>::DomainType &domain = levelSet->getDomain();

        std::vector<std::vector<T>> flagsReserve(
        levelSet->getNumberOfSegments());


        double pointsPerSegment =
        double(2 * levelSet->getDomain().getNumberOfPoints()) /
        double(levelSet->getLevelSetWidth());
        
        double time = 0.;

        double lsValSum = 0.;


        //shape
        if(method == 0){


        start = std::chrono::high_resolution_clock::now(); 

#pragma omp parallel num_threads((levelSet)->getNumberOfSegments())
        {
            int p = 0;
#ifdef _OPENMP
            p = omp_get_thread_num();
#endif


                curvaturShapeDerivatives1<T, D> curvatureCalculator(levelSet->getGrid().getGridDelta());

                //hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType> neighborStarIt(domain);


                auto &flagsSegment = flagsReserve[p];
                flagsSegment.reserve(pointsPerSegment);

                hrleVectorType<hrleIndexType, D> startVector =
                (p == 0) ? grid.getMinGridPoint()
                    : domain.getSegmentation()[p - 1];

                hrleVectorType<hrleIndexType, D> endVector =
                (p != static_cast<int>(domain.getNumberOfSegments() - 1))
                    ? domain.getSegmentation()[p]
                    : grid.incrementIndices(grid.getMaxGridPoint());


                for (hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType>
                    neighborIt(levelSet->getDomain(), startVector);
                    neighborIt.getIndices() < endVector; neighborIt.next()) {

                    if (!neighborIt.getCenter().isDefined() || std::abs(neighborIt.getCenter().getValue()) > 0.5) {
                        continue;
                    } 

                    //neighborStarIt.goToIndicesSequential(it.getStartIndices());

                    T curve = curvatureCalculator(neighborIt);

                    if(std::abs(curve) > flatBoundary){
                        flagsSegment.push_back(1);
                    }else{
                        flagsSegment.push_back(0);
                    }
                    
                }     
        }

        stop = std::chrono::high_resolution_clock::now(); 

        time = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();
        //general Formula
        }else if(method == 1){

        start = std::chrono::high_resolution_clock::now();

#pragma omp parallel num_threads((levelSet)->getNumberOfSegments())
            {
                int p = 0;
#ifdef _OPENMP
                p = omp_get_thread_num();
#endif

                curvaturGeneralFormula<T, D> curvatureCalculator(levelSet->getGrid().getGridDelta());

                //hrleCartesianPlaneIterator<hrleDomain<T, D>> planeIterator(levelSet->getDomain(), 1);

                auto &flagsSegment = flagsReserve[p];
                flagsSegment.reserve(pointsPerSegment);

                hrleVectorType<hrleIndexType, D> startVector =
                (p == 0) ? grid.getMinGridPoint()
                    : domain.getSegmentation()[p - 1];

                hrleVectorType<hrleIndexType, D> endVector =
                (p != static_cast<int>(domain.getNumberOfSegments() - 1))
                    ? domain.getSegmentation()[p]
                    : grid.incrementIndices(grid.getMaxGridPoint());

                //for (hrleSparseBoxIterator<typename lsDomain<T, D>::DomainType>
                //    neighborIt(levelSet->getDomain(), startVector, 1);
                //    neighborIt.getIndices() < endVector; neighborIt.next()) {

                for (hrleCartesianPlaneIterator<typename lsDomain<T, D>::DomainType>
                    neighborIt(levelSet->getDomain(), startVector, 1);
                    neighborIt.getIndices() < endVector; neighborIt.next()) {

                    if (!neighborIt.getCenter().isDefined() || std::abs(neighborIt.getCenter().getValue()) > 0.5) {
                        continue;
                    } 

                    T curve = curvatureCalculator(neighborIt);

                    //curveOutput.push_back(curve);

                    if(std::abs(curve) > flatBoundary){
                        flagsSegment.push_back(1);
                    }else{
                        flagsSegment.push_back(0);
                    }
                    
                } 

            }
            
            stop = std::chrono::high_resolution_clock::now(); 

            time = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();
        }

        for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i) 
            flaggedCells.insert(flaggedCells.end(),flagsReserve[i].begin(), flagsReserve[i].end());

        return time;
      

    }

};