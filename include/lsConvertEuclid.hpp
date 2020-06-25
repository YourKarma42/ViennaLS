#ifndef LS_CONVERT_EUCLID_HPP
#define LS_CONVERT_EUCLID_HPP

//TODO: paralellization

#include <lsPreCompileMacros.hpp>

#include <hrleSparseIterator.hpp>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>

#include <unordered_set>

//TODO: expand!
/// This class converts all level set values from a normalization with
/// Manhatten distance to a normalization with Euklidian distance.

template <class T, int D> class lsConvertEuclid {
  //typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  lsDomain<T, D> *levelSet = nullptr;

  std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> activePoints;

  //TODO: think of that value in euklidiant context 1/sqrt(2) ?
  T maxValue = 0.5;

public:
  lsConvertEuclid() {}

  lsConvertEuclid(lsDomain<T, D> &passedLevelSet, T passedMaxValue = 0.5)
      : levelSet(&passedLevelSet), maxValue(passedMaxValue) {
  }

  void setLevelSet(lsDomain<T, D> &passedLevelSet) {
    levelSet = &passedLevelSet;
  }

  void setMaxValue(const T passedMaxValue) { maxValue = passedMaxValue; }

  std::unordered_set<hrleVectorType<hrleIndexType, D>, typename hrleVectorType<hrleIndexType, D>::hash> &
  getActivePoints(){
      return activePoints;
  }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsConvertEuclid.")
          .print();
      return;
    }

    //increase lvl set size
    lsExpand<T, D>(*levelSet, 7).apply();


    double pointsPerSegment =
    double(2 * levelSet->getDomain().getNumberOfPoints()) /
    double(levelSet->getLevelSetWidth());

    auto grid = levelSet->getGrid();

    typename lsDomain<T, D>::DomainType &oldDomain = levelSet->getDomain();

    const T gridDelta = levelSet->getGrid().getGridDelta();

    //Initialize new level set with normalized values
    lsDomain<T,D> newLS(gridDelta);

    auto newGrid = newLS.getGrid();

    typename lsDomain<T, D>::DomainType &newDomain = newLS.getDomain();

    newDomain.initialize(oldDomain.getNewSegmentation(), oldDomain.getAllocation());

    //reserve space for active grid points
    activePoints.reserve(levelSet->getNumberOfPoints());


#pragma omp parallel num_threads(newDomain.getNumberOfSegments())
    {
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
       

    
        for (hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType>
            neighborIt(oldDomain, startVector);
            neighborIt.getIndices() < endVector; neighborIt.next()) {

            auto &centerIt = neighborIt.getCenter();
            if (!centerIt.isDefined() || std::abs(centerIt.getValue()) > 0.5) {
            //write undefined run in new level set
            newDomain.getDomainSegment(p).insertNextUndefinedPoint(neighborIt.getIndices(), 
            (centerIt.getValue()<0.) ? lsDomain<T, D>::NEG_VALUE : lsDomain<T, D>::POS_VALUE);

            continue;
            } 

            //TODO: probably remove
            //TODO: maby there is an error here when paralell check calculate normal vectors
            activePoints.insert(centerIt.getStartIndices());


            //TODO: debug save old LS value for comparison later
            T oldLSValue = centerIt.getValue();

            //calculate normal vector of defined grid point
            std::array<T, D> n;
            T normN = 0.;
            for (int i = 0; i < D; i++) {
            
            T pos1 = neighborIt.getNeighbor(i).getValue();
            T neg1 = neighborIt.getNeighbor(i + D).getValue();
            n[i] = (pos1 - neg1) * 0.5;
            
            normN += n[i] * n[i];
            }
            normN = std::sqrt(normN);

            for (int i = 0; i < D; i++) {
            n[i] /= normN;
            }

            //renormalized ls value in direction of the normal vector
            T max = 0.;
            for (unsigned i = 0; i < D; ++i) {

                if (std::abs(n[i]) > max) {
                    max = std::abs(n[i]);
                }
            }
    
            T newLsValue = centerIt.getValue() * max;

            newDomain.getDomainSegment(p).insertNextDefinedPoint(neighborIt.getIndices(), newLsValue);

        }
    }

    newDomain.finalize();

  }


};







#endif // LS_CONVERT_EUCLID_HPP