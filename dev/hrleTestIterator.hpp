#ifndef HRLE_TEST_ITERATOR_HPP
#define HRLE_TEST_ITERATOR_HPP

//TODO: change location into HRLE
#include "../../ViennaHRLE/include/hrleSparseIterator.hpp"
#include "../../ViennaHRLE/include/hrleSparseOffsetIterator.hpp"


template <class hrleDomain> class hrleTestIterator {




    static constexpr int D = hrleDomain::dimension;

    hrleDomain &domain;
    const unsigned size;
    const hrleIndexType centerIndex;
    hrleVectorType<hrleIndexType, D> currentCoords;
    std::vector<hrleSparseOffsetIterator<hrleDomain>> neighborIterators;



    bool isDefined() const {
        if (neighborIterators[centerIndex].isDefined())
        return true;
        for (int i = 0; i < 2 * D; i++) {
        if (neighborIterators[i].isDefined())
            return true;
        }
        return false;
    }

    bool isCenterDefined() const {
        if (neighborIterators[centerIndex].isDefined())
            return true;

        return false;
    }



    /// push offset iterators lexicographically into std::vector from -order to
    /// +order
    template <class V> void initializeNeigbors(const V &v) {

        //TODO: rewrite!!!
        const unsigned numNeighbors = size;

        neighborIterators.reserve(numNeighbors);

        for (unsigned i = 0; i < numNeighbors; ++i) {
        hrleVectorType<hrleIndexType, D> offset = indexToCoordinate(i);
        neighborIterators.push_back(
            hrleSparseOffsetIterator<hrleDomain>(domain, offset, v));
        }
    }

    //TODO: rerwirte currently only for 3D only allocates fixed stencil
    hrleVectorType<hrleIndexType, D>
    indexToCoordinate(hrleIndexType index) const {
        hrleVectorType<hrleIndexType, D> coordinate;

        unsigned order = 1;

        coordinate[1] = index / 3;
        coordinate[0] = index % 3;

        // shift to the middle
        for (unsigned i = 0; i < D; ++i)
        coordinate[i] -= order;

        return coordinate;
    }

    //TODO: rerwirte currently only for 3D only allocates fixed stencil
    template <class V> hrleIndexType coordinateToIndex(V coordinate) const {
        unsigned order = 1;
        // shift to the middle
        for (unsigned i = 0; i < D; ++i)
        coordinate[i] += order;

        hrleIndexType index = 0;
        if (D > 2) {
            index += coordinate[2] * 9;
        }
        index += coordinate[1] * 3;
        index += coordinate[0];

        return index;
    }

    // make post in/decrement private, since they should not be used, due to the
    // size of the structure
    hrleTestIterator<hrleDomain>
    operator++(int) { // use pre increment instead
        return *this;
    }
    hrleTestIterator<hrleDomain>
    operator--(int) { // use pre decrement instead
        return *this;
    }

    public:

    hrleTestIterator(hrleDomain &passedDomain,
                            const hrleVectorType<hrleIndexType, D> &v,
                            const unsigned size = 1)
        : domain(passedDomain), size(size),
            centerIndex(coordinateToIndex(hrleVectorType<hrleIndexType, D>(0))),
            currentCoords(v) {

        initializeNeigbors(v);
    }

    hrleTestIterator(hrleDomain &passedDomain,
                            const unsigned size = 1)
        : domain(passedDomain), size(size),
            centerIndex(coordinateToIndex(hrleVectorType<hrleIndexType, D>(0))),
            currentCoords(domain.getGrid().getMinGridPoint()) {

        initializeNeigbors(passedDomain.getGrid().getMinIndex());
    }
/*
    hrleTestIterator<hrleDomain> &operator++() {
        next();
        return *this;
    }

    hrleTestIterator<hrleDomain> &operator--() {
        previous();
        return *this;
    }
*/

    hrleSparseOffsetIterator<hrleDomain> &getCenter() {
        return neighborIterators[centerIndex];
    }

    const hrleVectorType<hrleIndexType, D> &getIndices() { return currentCoords; }

    unsigned getSize() { return neighborIterators.size(); }

    bool isFinished() const { return getCenter().isFinished(); }

    //TODO:next and previous are not implemented yet

    template <class V> void goToIndices(V &v) {
        const unsigned numNeighbors = neighborIterators.size();
        getCenter().goToIndices(v);
        for (int i = 0; i < 1; ++i) {
            for (int j = 0; j < numNeighbors; ++j) {
                neighborIterators[j].goToIndices(v);
            }
        }
    }


};



#endif // HRLE_TEST_ITERATOR_HPP