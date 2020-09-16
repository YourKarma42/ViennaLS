#ifndef HRLE_TEST_ITERATOR_HPP
#define HRLE_TEST_ITERATOR_HPP

//TODO: change location into HRLE
//#include "../../../ViennaHRLE/viennahrle/include/hrleSparseIterator.hpp"
//#include "../../../ViennaHRLE/viennahrle/include/hrleSparseOffsetIterator.hpp"

#include "../../ViennaHRLE/include/hrleSparseIterator.hpp"
#include "../../ViennaHRLE/include/hrleSparseOffsetIterator.hpp"

/// This neighbor iterator consists of (2*order+1)^dimension
/// hrleSparseOffsetIterator s for the cartesian neighbors and the center.
/// Whenever one of these iterators reach a defined grid point, the square
/// iterator stops.
/// Neighbors are indexed lexicographically from negative cartesian directions:
/// order 1:            order 2:d
///                     20 21 22 23 24
/// 6 7 8               15 16 17 18 19
/// 3 4 5               10 11 12 13 14
/// 0 1 2               5  6  7  8  9
///                     0  1  2  3  4
/// center: 4           center: 12
template <class hrleDomain> class hrleCartesianPlaneIterator {

  typedef typename std::conditional<std::is_const<hrleDomain>::value,
                                    const typename hrleDomain::hrleValueType,
                                    typename hrleDomain::hrleValueType>::type
      hrleValueType;

  static constexpr int D = hrleDomain::dimension;

  hrleDomain &domain;
  const unsigned order;
  const hrleIndexType sideLength;
  const hrleIndexType sliceArea;
  const hrleIndexType centerIndex;
  hrleSparseIterator<hrleDomain> centerIterator;

  hrleVectorType<hrleIndexType, D> currentCoords;
  std::vector<hrleSparseOffsetIterator<hrleDomain>> neighborIterators;

  bool isDefined() const {
    if (centerIterator.isDefined())
      return true;
    //TODO: rethink!
    for (int i = 0; i <  D; i++) {
      if (neighborIterators[i].isDefined())
        return true;
    }
    return false;
  }

  unsigned axisEnd = (D * 2 * order);
  unsigned planeEnd = 4 * order * order;

  // | Star Iterator |  YZ plane  |  XZ plane  |  XY plane  |
  // | D * 2 * order |  4*order^2 |  4*order^2 |  4*order^2 |

  hrleVectorType<hrleIndexType, D>
  indexToCoordinate(hrleIndexType index) const {

    hrleVectorType<hrleIndexType, D> coordinate(hrleIndexType(0));

    //TESTING
    //std::cout << index << std::endl;

    if(index < axisEnd){
      //The index is part of the star iterator    

      unsigned plane = index / (2*order);

      unsigned axisPos = index % (2*order);

      if(axisPos/order == 0){
        coordinate[plane] = -(axisPos/order + order); 
      }else{
        //TODO:only order 1
        coordinate[plane] = (axisPos/order); 
      }

    }else{
      //shift coordinates to planes

      index -= axisEnd;

      unsigned zeroCoordinate = index / planeEnd;

      index = index % planeEnd;

      unsigned depthAxis = index / (2*order);

      unsigned sideAxis = index % (2*order);


      unsigned pos1 = (zeroCoordinate + 1) % D;
      unsigned pos2 = (zeroCoordinate + 2) % D;
      //flip axis for the last 2 indices so that we always look at the plane from the lower dimesnion
      if((zeroCoordinate + 2) % D == 0){
        unsigned tmp = pos1;
        pos1 = pos2;
        pos2 = tmp;
      }


      if(depthAxis >= order){
        coordinate[pos1] = order; 
      }else{
        coordinate[pos1] = -order; 
      }

      if(sideAxis >= order){
        coordinate[pos2] = order; 
      }else{
        coordinate[pos2] = -order; 
      }

    }

    //TESTING
    //std::cout << coordinate << "  " << coordinateToIndex(coordinate) << std::endl;
 
    return coordinate;
  }

  template <class V> hrleIndexType coordinateToIndex(V coordinate) const {

    //TODO:consider (0,0,0) vector

    std::vector<int> pos(0);
    int anz = 0;

    unsigned index = 0;

    unsigned zeroPos = 0;

    for(int i = 0; i<D; i++){

      if(coordinate[i] != 0){      
        pos.push_back(i);
        anz++;
      }else{
        zeroPos = i;
      }

    }

    if(anz == 1){
      if(coordinate[pos[0]] > 0){
        return (pos[0]*2*order)+1;
      }else{
        return (pos[0]*2*order);
      }
    }else{

      index = axisEnd + (zeroPos * planeEnd);

      if(coordinate[pos[0]] > 0){
        index += 2*order * (order + coordinate[pos[0]] - 1);
      }else{
        index += 2*order * (order + coordinate[pos[0]]);
      }

      if(order == 1){
        if(coordinate[pos[1]] > 0){
          index += coordinate[pos[1]];
        }else{
          index += order + coordinate[pos[1]];
        }
      }else{
        //TODO test for order = 2

      }

    }


    return index;
  }

  /// push offset iterators lexicographically into std::vector from -order to
  /// +order
  template <class V> void initializeNeigbors(const V &v) {

    //TODO for D = 2 this is equivalent to a box iterator star iterator!
    const unsigned numNeighbors = unsigned(std::pow((1 + 2 * order), D)-(8*std::pow(order,D)));

    neighborIterators.reserve(numNeighbors - 1);

    for (unsigned i = 0; i < numNeighbors - 1; ++i) {
      hrleVectorType<hrleIndexType, D> offset = indexToCoordinate(i);
      neighborIterators.push_back(
          hrleSparseOffsetIterator<hrleDomain>(domain, offset, v));
    }
  }

  // make post in/decrement private, since they should not be used, due to the
  // size of the structure
  hrleCartesianPlaneIterator<hrleDomain>
  operator++(int) { // use pre increment instead
    return *this;
  }
  hrleCartesianPlaneIterator<hrleDomain>
  operator--(int) { // use pre decrement instead
    return *this;
  }

public:
  hrleCartesianPlaneIterator(hrleDomain &passedDomain,
                        const hrleVectorType<hrleIndexType, D> &v,
                        const unsigned passedOrder = 1)
      : domain(passedDomain), order(passedOrder),
        sideLength(1 + 2 * passedOrder), sliceArea(sideLength * sideLength),
        centerIndex(1),
        currentCoords(v),
        centerIterator(passedDomain, v) {

    initializeNeigbors(v);
  }

  hrleCartesianPlaneIterator(hrleDomain &passedDomain,
                        const unsigned passedOrder = 1)
      : domain(passedDomain), order(passedOrder),
        sideLength(1 + 2 * passedOrder), sliceArea(sideLength * sideLength),
        centerIndex(1),
        currentCoords(domain.getGrid().getMinGridPoint()),
        centerIterator(passedDomain) {

    initializeNeigbors(passedDomain.getGrid().getMinIndex());
  }

  hrleCartesianPlaneIterator<hrleDomain> &operator++() {
    next();
    return *this;
  }

  hrleCartesianPlaneIterator<hrleDomain> &operator--() {
    previous();
    return *this;
  }

  void next() {
    const int numNeighbors = int(neighborIterators.size());
    std::vector<bool> increment(numNeighbors + 1, false);
    increment[numNeighbors] = true;

    hrleVectorType<hrleIndexType, D> end_coords =
        neighborIterators[centerIndex].getEndIndices();
    for (int i = 0; i < numNeighbors; i++) {
      if (i == centerIndex)
        continue;

      switch (compare(end_coords, neighborIterators[i].getEndIndices())) {
      case 1:
        end_coords = neighborIterators[i].getEndIndices();
        increment = std::vector<bool>(numNeighbors + 1, false);
      case 0:
        increment[i] = true;
      }
    }

    if (increment[numNeighbors])
      neighborIterators[centerIndex].next();
    for (int i = 0; i < numNeighbors; i++)
      if (increment[i])
        neighborIterators[i].next();

    currentCoords = domain.getGrid().incrementIndices(end_coords);
  }

  void previous() {
    const int numNeighbors = neighborIterators.size();
    std::vector<bool> decrement(numNeighbors + 1, false);
    decrement[numNeighbors] = true;

    hrleVectorType<hrleIndexType, D> start_coords =
        neighborIterators[centerIndex].getStartIndices();
    for (int i = 0; i < numNeighbors; i++) {
      if (i == centerIndex)
        continue;
      switch (compare(start_coords, neighborIterators[i].getStartIndices())) {
      case -1:
        start_coords = neighborIterators[i].getStartIndices();
        decrement = std::vector<bool>(numNeighbors + 1, false);
      case 0:
        decrement[i] = true;
      }
    }

    if (decrement[numNeighbors])
      neighborIterators[centerIndex].previous();
    for (int i = 0; i < numNeighbors; i++)
      if (decrement[i])
        neighborIterators[i].previous();

    currentCoords = domain.getGrid().decrementIndices(start_coords);
  }

  hrleSparseOffsetIterator<hrleDomain> &getNeighbor(int index) {
    return neighborIterators[index];
  }

  hrleSparseOffsetIterator<hrleDomain> &getNeighbor(unsigned index) {
    return neighborIterators[index];
  }

  template <class V>
  hrleSparseOffsetIterator<hrleDomain> &getNeighbor(V relativeCoordinate) {
    return neighborIterators[coordinateToIndex(relativeCoordinate)];
  }

hrleSparseIterator<hrleDomain> &getCenter() { return centerIterator; }

  const hrleVectorType<hrleIndexType, D> &getIndices() { return currentCoords; }

  unsigned getSize() { return neighborIterators.size(); }

  bool isFinished() const { return getCenter().isFinished(); }

  /// Sets the iterator to position v.
  /// Uses random access to move, so it is be slower
  /// than goToIndicesSequential for repeated serial calls.
  template <class V> void goToIndices(V &v) {
    const unsigned numNeighbors = neighborIterators.size();
    getCenter().goToIndices(v);
    //for (int i = 0; i < int(order); ++i) {
    for (int j = 0; j < numNeighbors; ++j) {
      neighborIterators[j].goToIndices(v);
    }
    //}
  }

  /// Advances the iterator to position v.
  /// If v is lexicographically higher than the current position
  /// the iterator will be moved back to v.
  /// If v is lexicographically smaller than the current position
  /// then the iterator will be moved until it reaches v
  template <class V> void goToIndicesSequential(const V &v) {
    if (v >= currentCoords) {
      while (v > currentCoords) {
        next();
      }
    } else {
      while (v < currentCoords) {
        previous();
      }
    }
  }
};

template <class hrleDomain>
using hrleConstCartesianPlaneIterator = hrleCartesianPlaneIterator<const hrleDomain>;

#endif // HRLE_TEST_ITERATOR_HPP