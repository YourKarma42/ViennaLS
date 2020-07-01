#include <iostream>

#include <lsDomain.hpp>

#include <hrleSparseBoxIterator.hpp>



template <class T, int D> class curvaturShapeDerivatives1{

    private:

    //hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator;

    T gridDelta;

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

    curvaturShapeDerivatives1(T mGD)
    :  gridDelta(mGD){
    }
    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */

    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){

        //calculate all needed derivatives in the xy yz and xz plane

        std::array<T, 3> centralDiff;

        std::array<T, 6> oneSidedeDiff;

        std::array<T, 6> derivativesPos;

        std::array<T, 6> derivativesNeg;
        

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            hrleVectorType<hrleIndexType, D> test = neighborIterator.getIndices();

            //TODO: ask people if there is a more elegant solution
            int first_pos  = i;
            int second_pos = (i+1) % D;
            

            posUnit[first_pos] = 1;
            negUnit[first_pos] = -1;
        
            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            //phi_0 = convertToFuncVal(phi_0, neighborIterator.getCenter().getStartIndices(), gridDelta);

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            //phi_px = convertToFuncVal(phi_px, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
            //phi_nx = convertToFuncVal(phi_nx, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

            //This is necessary to to keep grid axis directions consistent in each slice
            if(i > 0){
                // probably not evrything needed
                posUnit[first_pos] = 0;
                negUnit[first_pos] = 0;

                first_pos  = second_pos;
                second_pos = i;

                posUnit[first_pos] = 1;
                negUnit[first_pos] = -1;
            }

            posUnit[second_pos] = 1;
            negUnit[second_pos] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            //phi_pp = convertToFuncVal(phi_pp, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
            //phi_np = convertToFuncVal(phi_np, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

            posUnit[second_pos] = -1;
            negUnit[second_pos] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            //phi_pn = convertToFuncVal(phi_pn, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
            //phi_nn = convertToFuncVal(phi_nn, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

            posUnit[first_pos] = 0;
            negUnit[first_pos] = 0;

            posUnit[second_pos] = 1;
            negUnit[second_pos] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

            //phi_py = convertToFuncVal(phi_py, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
            //phi_ny = convertToFuncVal(phi_ny, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

            //central
            centralDiff[i] = (phi_px - phi_nx)*0.5;

            //one sided
            oneSidedeDiff[i] = phi_px - phi_0;
            oneSidedeDiff[i+3] = phi_nx - phi_0;

            //central outer
            derivativesPos[i]   = (phi_pp - phi_np)*0.5; 
            derivativesPos[i+3] = (phi_pp - phi_pn)*0.5; 

            derivativesNeg[i] = (phi_pn- phi_nn)*0.5; 
            derivativesNeg[i+3] = (phi_np - phi_nn)*0.5; 

        

    }


};


template <class T, int D> class curvatur1{

    private:

    //hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator;

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

    curvatur1(T mGD)
    :  gridDelta(mGD){
    }

    //curvatur1(hrleSparseBoxIterator<hrleDomain<T, D>> &ni, T mGD)
    //: neighborIterator(ni), gridDelta(mGD){
    //}

    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */

    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){

        //calculate all needed derivatives in the xy yz and xz plane

        std::array<T, 3> centralDiff;

        std::array<T, 6> oneSidedeDiff;

        std::array<T, 6> derivativesPos;

        std::array<T, 6> derivativesNeg;
        

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            hrleVectorType<hrleIndexType, D> test = neighborIterator.getIndices();

            //TODO: ask people if there is a more elegant solution
            int first_pos  = i;
            int second_pos = (i+1) % D;
            

            posUnit[first_pos] = 1;
            negUnit[first_pos] = -1;
        
            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            //phi_0 = convertToFuncVal(phi_0, neighborIterator.getCenter().getStartIndices(), gridDelta);

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            //phi_px = convertToFuncVal(phi_px, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
            //phi_nx = convertToFuncVal(phi_nx, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

            //This is necessary to to keep grid axis directions consistent in each slice
            if(i > 0){
                // probably not evrything needed
                posUnit[first_pos] = 0;
                negUnit[first_pos] = 0;

                first_pos  = second_pos;
                second_pos = i;

                posUnit[first_pos] = 1;
                negUnit[first_pos] = -1;
            }

            posUnit[second_pos] = 1;
            negUnit[second_pos] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            //phi_pp = convertToFuncVal(phi_pp, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
            //phi_np = convertToFuncVal(phi_np, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

            posUnit[second_pos] = -1;
            negUnit[second_pos] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            //phi_pn = convertToFuncVal(phi_pn, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
            //phi_nn = convertToFuncVal(phi_nn, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

            posUnit[first_pos] = 0;
            negUnit[first_pos] = 0;

            posUnit[second_pos] = 1;
            negUnit[second_pos] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

            //phi_py = convertToFuncVal(phi_py, neighborIterator.getNeighbor(posUnit).getStartIndices(), gridDelta);
            //phi_ny = convertToFuncVal(phi_ny, neighborIterator.getNeighbor(negUnit).getStartIndices(), gridDelta);

            //central
            centralDiff[i] = (phi_px - phi_nx)*0.5;

            //one sided
            oneSidedeDiff[i] = phi_px - phi_0;
            oneSidedeDiff[i+3] = phi_nx - phi_0;

            //central outer
            derivativesPos[i]   = (phi_pp - phi_np)*0.5; 
            derivativesPos[i+3] = (phi_pp - phi_pn)*0.5; 

            derivativesNeg[i] = (phi_pn- phi_nn)*0.5; 
            derivativesNeg[i+3] = (phi_np - phi_nn)*0.5; 

        }

        
            //n_x_pos - n_x_neg
        T n_x = (oneSidedeDiff[0] /
                (2*(std::sqrt(oneSidedeDiff[0]*oneSidedeDiff[0]) + 
                (derivativesPos[3] + centralDiff[1]*0.5 + 
                (derivativesPos[5] + centralDiff[2])*0.5)))) 
                -
                (oneSidedeDiff[3] /
                (2*(std::sqrt(oneSidedeDiff[3]*oneSidedeDiff[3]) + 
                (derivativesNeg[3] + centralDiff[1]*0.5 + 
                (derivativesNeg[5] + centralDiff[2])*0.5))));

        T n_y = (oneSidedeDiff[1] /
                (2*(std::sqrt(oneSidedeDiff[1]*oneSidedeDiff[1]) + 
                (derivativesPos[0] + centralDiff[0]*0.5 + 
                (derivativesPos[1] + centralDiff[2])*0.5)))) 
                -
                (oneSidedeDiff[4] /
                (2*(std::sqrt(oneSidedeDiff[4]*oneSidedeDiff[4]) + 
                (derivativesNeg[0] + centralDiff[0]*0.5 + 
                (derivativesNeg[1] + centralDiff[2])*0.5))));

        T n_z = (oneSidedeDiff[2] /
                (2*(std::sqrt(oneSidedeDiff[2]*oneSidedeDiff[2]) + 
                (derivativesPos[2] + centralDiff[0]*0.5 + 
                (derivativesPos[4] + centralDiff[1])*0.5)))) 
                -
                (oneSidedeDiff[5] /
                (2*(std::sqrt(oneSidedeDiff[5]*oneSidedeDiff[5]) + 
                (derivativesNeg[2] + centralDiff[0]*0.5 + 
                (derivativesNeg[4] + centralDiff[1])*0.5))));
                
            
        return n_x + n_y + n_z;
        

    }


};
