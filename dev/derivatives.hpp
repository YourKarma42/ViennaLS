#include <iostream>

#include <lsDomain.hpp>

#include <hrleSparseBoxIterator.hpp>




void test1(){
    std::cout << "blubdiblala" << std::endl;
}

/*
template <class T, int D> class curveCaluclator {



    public:

    virtual void apply() = 0;




};*/




template <class T, int D> class curvatur1{

    private:

    hrleSparseBoxIterator<hrleDomain<T, D>> neighborIterator;

    void calcDerivatives(){

    }

    public:

    curvatur1(hrleSparseBoxIterator<hrleDomain<T, D>> &ni)
    : neighborIterator(ni){
    }

    void apply(){

    }

};
