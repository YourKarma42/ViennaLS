#include <iostream>

#include <lsDomain.hpp>

#include <hrleSparseBoxIterator.hpp>



template <class T, int D> class silvacoLikeFlagging{

    private:

    lsDomain<T, D> *levelSet = nullptr;

    public:

    silvacoLikeFlagging(lsDomain<T, D> &passedLevelSet)
      : levelSet(&passedLevelSet){
    }


    void apply(){

        //calculate normals of the lvl set

        //Iterate over level set and flag cells

       
   

    }

};