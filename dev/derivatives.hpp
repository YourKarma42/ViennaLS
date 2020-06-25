#include <iostream>

#include <lsDomain.hpp>


typedef double NumericType;


void test1(){
    std::cout << "blubdiblala" << std::endl;
}

class curveCaluclator {


    private:

    hrleSparseBoxIterator<hrleDomain<NumericType, D>> neighborIterator;

    public:

    curveCaluclator(hrleSparseBoxIterator<hrleDomain<NumericType, D>> ni):
    neighborIterator(ni){

    }

    void apply(){

    }




}