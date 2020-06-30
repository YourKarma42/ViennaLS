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

    hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator;

    T gridDelta;

    void calcDerivatives(){
        

    }

    T convertToFuncVal(T lsVal, hrleVectorType<hrleIndexType, D> gridPoint, T gridDelta){

        T distGridPoint = 0.;

        for(int i = 0; i < D; i++){
            distGridPoint += (gridPoint[i]*gridDelta) * (gridPoint[i]*gridDelta);
        }

        distGridPoint = std::sqrt(distGridPoint);

        if(lsVal > 0.){
            return (lsVal*gridDelta) + distGridPoint;
        }else{
            return ((lsVal*gridDelta) + distGridPoint);
        }

    }

    public:

    curvatur1(hrleSparseBoxIterator<hrleDomain<T, D>> &ni, T mGD)
    : neighborIterator(ni), gridDelta(mGD){
    }

    void apply(){

        //calculate all needed derivatives in the xy yz and xz plane

        std::array<T, 3> centralDiff;

        std::array<T, 6> oneSidedeDiff;

        std::array<T, 12> derivatives;
        

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            hrleVectorType<hrleIndexType, D> test = neighborIterator.getIndices();

            posUnit[i] = 1;
            negUnit[i] = -1;

            int second_pos = (i+1) % D;

            if(i == 1){
                second_pos = 0;
            }

            //get required ls values
            NumericType phi_0 = neighborIterator.getCenter().getValue();

            phi_0 = convertToFuncVal(phi_0, neighborIterator.getCenter().getStartIndices(), gridDelta);

            NumericType phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            NumericType phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            phi_px = convertToFuncVal(phi_px, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
            phi_nx = convertToFuncVal(phi_nx, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

            posUnit[second_pos] = 1;
            negUnit[second_pos] = 1;

            NumericType phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            NumericType phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            phi_pp = convertToFuncVal(phi_pp, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
            phi_np = convertToFuncVal(phi_np, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

            posUnit[second_pos] = -1;
            negUnit[second_pos] = -1;

            NumericType phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            NumericType phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            phi_pn = convertToFuncVal(phi_pn, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
            phi_nn = convertToFuncVal(phi_nn, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

            posUnit[i] = 0;
            negUnit[i] = 0;

            posUnit[second_pos] = 1;
            negUnit[second_pos] = -1;

            NumericType phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            NumericType phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

            phi_py = convertToFuncVal(phi_py, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
            phi_ny = convertToFuncVal(phi_ny, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

            //central
            centralDiff[i] = (phi_px - phi_nx)*0.5;

            //one sided
            oneSidedeDiff[i + 1] = phi_px - phi_0;
            oneSidedeDiff[i + 2] = phi_nx - phi_0;

            //central outer
            derivatives[i*4]     = (phi_pp - phi_np)*0.5; //D_i^{i+1%3}
            derivatives[i*4 + 1] = (phi_pn - phi_nn)*0.5; //D_i^{i-1%3}
            derivatives[i*4 + 2] = (phi_pp - phi_pn)*0.5; //D_i^{+i+1%3}
            derivatives[i*4 + 3] = (phi_np - phi_nn)*0.5; //D_i^{+i+1%3}

        }

        for(int i = 0; i < D; i++) {

            T n_xP = oneSidedeDiff[i] /
                (2*(std::sqrt(oneSidedeDiff[i]*oneSidedeDiff[i]) + 
                (derivatives[1] + centralDiff[(i+1)%D])*0.5 + 
                (derivatives[1] + centralDiff[(i+2)%D])*0.5));
            

        }

    }


};
