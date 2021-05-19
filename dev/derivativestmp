#include <iostream>

#include <lsDomain.hpp>

#include <hrleSparseBoxIterator.hpp>

#include "hrleTestIterator.hpp"

//TODO: calculate grid delta values in consturctor

template <class T, int D> class baseDerivative{
    public:
    virtual T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator) {return -1;}
    virtual T operator()(hrleCartesianPlaneIterator<hrleDomain<T, D>> & neighborIterator) {return -1;}
};

template <class T, int D> class curvaturGeneralFormula : public baseDerivative<T, D>{

    private:

    //hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator;

    T gridDelta;

    T twoGD = 0;
    T gDSq = 0;
    T fourGDsq = 0;

    public:

    curvaturGeneralFormula(T mGD)
    :  gridDelta(mGD){
        twoGD = 1./(2.*gridDelta);
        gDSq = 1./(gridDelta*gridDelta);
        fourGDsq = 1./(4.*gridDelta*gridDelta);


    }
    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */
    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator) {

        //calculate all needed derivatives in the xy yz and xz plane



        std::array<T, 9> d;
       
        //get required ls values

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            posUnit[i] = 1;
            negUnit[i] = -1;

            int second_pos = (i+1) % D;

            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_pos] = 1;
            negUnit[second_pos] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_pos] = -1;
            negUnit[second_pos] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[i] = 0;
            negUnit[i] = 0;

            posUnit[second_pos] = 1;
            negUnit[second_pos] = -1;

            // first order derivative
            d[i] = (phi_px - phi_nx)*twoGD;

            d[i+3] = (phi_px - 2.*phi_0 + phi_nx)*gDSq;

            d[i+6] = (phi_pp - phi_pn -phi_np + phi_nn)*fourGDsq;

        }

        T norm_grad_pow3 = 0.;
        for(int i = 0; i < D; i++){

            norm_grad_pow3 += d[i]*d[i];
        }
        norm_grad_pow3 = std::sqrt(norm_grad_pow3);

        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;

        if(D == 2){
            return (d[3]*d[1]*d[1] - 2.*d[1]*d[0]*d[6] + d[4]*d[0]*d[0])
                    /(norm_grad_pow3);
        }


        //expanded fom of the equation
        return
        //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx + F_yy)
        (d[0]*d[0]*(d[4] + d[5]) + d[1]*d[1]*(d[3] + d[5]) + d[2]*d[2]*(d[3] + d[4]) +

        //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
        -2.*(d[0]*d[1]*d[6] + d[0]*d[2]*d[8] + d[1]*d[2]*d[7]))
                
        // /2*(F_x² + F_y² + F_z²)^(3/2)
        /(2.*norm_grad_pow3);

        //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)


    }

    //T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){
    T operator()(hrleCartesianPlaneIterator<hrleDomain<T, D>> & neighborIterator) {

        //calculate all needed derivatives in the xy yz and xz plane

        std::array<T, 9> d;
       
        //get required ls values

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            posUnit[i] = 1;
            negUnit[i] = -1;

            int second_pos = (i+1) % D;

            //if(i == 1){
            //    second_pos = 0;
            //}

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_pos] = 1;
            negUnit[second_pos] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_pos] = -1;
            negUnit[second_pos] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[i] = 0;
            negUnit[i] = 0;

            posUnit[second_pos] = 1;
            negUnit[second_pos] = -1;

            // first order derivative
            d[i] = (phi_px - phi_nx)*twoGD;

            d[i+3] = (phi_px - 2.*phi_0 + phi_nx)*gDSq;

            d[i+6] = (phi_pp - phi_pn -phi_np + phi_nn)*fourGDsq;

        }


        T norm_grad_pow3 = 0.;
        for(int i = 0; i < D; i++){

            norm_grad_pow3 += d[i]*d[i];
        }
        norm_grad_pow3 = std::sqrt(norm_grad_pow3);

        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;

        if(D == 2){
            return (d[3]*d[1]*d[1] - 2.*d[1]*d[0]*d[6] + d[4]*d[0]*d[0])
                    /(norm_grad_pow3);
        }

        //expanded fom of the equation
        return
        //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx + F_yy)
        (d[0]*d[0]*(d[4] + d[5]) + d[1]*d[1]*(d[3] + d[5]) + d[2]*d[2]*(d[3] + d[4]) +

        //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
        -2.*(d[0]*d[1]*d[6] + d[0]*d[2]*d[8] + d[1]*d[2]*d[7]))
                
        // /2*(F_x² + F_y² + F_z²)^(3/2)
        /(2.*norm_grad_pow3);

        //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)


    }

};


template <class T, int D> class curvaturShapeDerivatives1 : public baseDerivative<T, D>{

    private:

    T gridDelta;

    T GDSQ = 0;

    public:

    curvaturShapeDerivatives1(T mGD)
    :  gridDelta(mGD){
        GDSQ = 1./(gridDelta*gridDelta);
    }
    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */
    T operator()(hrleCartesianPlaneIterator<hrleDomain<T, D>> & neighborIterator) {return -1;}

    T operator()(hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType> & neighborIterator){

        std::array<T, 3> scondOrderDerivatives;
       
        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);
        
            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(i).getValue();
            T phi_nx = neighborIterator.getNeighbor(i+D).getValue();

            scondOrderDerivatives[i] = (phi_px - 2.*phi_0 + phi_nx)*GDSQ;  
        }

        T result = 0.;

        for(int i=0; i<D ; i++)
            result += scondOrderDerivatives[i];
        //mean curvature is trace devided by 2
        return result*0.5;

    }


    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){

        //calculate all needed derivatives in the xy yz and xz plane

        std::array<T, 3> scondOrderDerivatives;
       

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            //TODO: ask people if there is a more elegant solution
            int first_pos  = i;
            //int second_pos = (i+1) % D;
            

            posUnit[first_pos] = 1;
            negUnit[first_pos] = -1;
        
            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();


            //central
            scondOrderDerivatives[i] = (phi_px - 2.*phi_0 + phi_nx)*GDSQ;  

        }


        T result = 0.;

        for(int i=0; i<D ; i++)
            result += scondOrderDerivatives[i];
        //mean curvature is trace devided by 2
        return result*0.5;

    }

};



template <class T, int D> class curvaturShapeDerivatives2 : public baseDerivative<T, D>{

    private:

    T gridDelta;

    T threeGDsq = 0;

    public:

    curvaturShapeDerivatives2(T mGD)
    :  gridDelta(mGD){
        threeGDsq = 1./(3.*gridDelta*gridDelta);
    }
    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */

   

    T operator()(hrleCartesianPlaneIterator<hrleDomain<T, D>> & neighborIterator){

        //calculate all needed derivatives in the xy yz and xz plane

        std::array<T, D> scondOrderDerivatives;

        for(int i = 0; i < D; i ++){
       
            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_pos = i;
            int second_pos = (i+1)%D;
            

            posUnit[first_pos] = 1;
            negUnit[first_pos] = -1;
        
            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_pos] = 1;
            negUnit[second_pos] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();


            posUnit[second_pos] = -1;
            negUnit[second_pos] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_pos] = 0;
            negUnit[first_pos] = 0;

            posUnit[second_pos] = 1;
            negUnit[second_pos] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();


            //central
            scondOrderDerivatives[i] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*threeGDsq;
        }

        T result = 0.;

        for(int i=0; i<D ; i++)
            result += scondOrderDerivatives[i];
        //mean curvature is trace devided by 2
        return result*0.5;

    }

    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){

        //calculate all needed derivatives in the xy yz and xz plane

        std::array<T, D> scondOrderDerivatives;

        for(int i = 0; i < D; i ++){
       
            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_pos = i;
            int second_pos = (i+1)%D;
            

            posUnit[first_pos] = 1;
            negUnit[first_pos] = -1;
        
            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_pos] = 1;
            negUnit[second_pos] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();


            posUnit[second_pos] = -1;
            negUnit[second_pos] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_pos] = 0;
            negUnit[first_pos] = 0;

            posUnit[second_pos] = 1;
            negUnit[second_pos] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

            scondOrderDerivatives[i] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*threeGDsq;
        }

        T result = 0.;

        for(int i=0; i<D ; i++)
            result += scondOrderDerivatives[i];
        //mean curvature is trace devided by 2
        return result*0.5;

    }

};

template <class T, int D> class curvaturShapeBias : public baseDerivative<T, D>{

    private:

    T gridDelta;

    T twoGD = 0;
    T threeGDsq = 0;

    public:

    curvaturShapeBias(T mGD)
    :  gridDelta(mGD){
        twoGD = 1./(2.*gridDelta);
        threeGDsq = 1./(3.*gridDelta*gridDelta);
    }
    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */

   

    T operator()(hrleCartesianPlaneIterator<hrleDomain<T, D>> & neighborIterator){


        //calculate Gradient

        std::array<T,D> gradient;
        //TODO: think of using upwind normals?
        for(int i = 0; i < D; i++){

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            posUnit[i] = 1;
            negUnit[i] = -1;

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            gradient[i] = std::abs((phi_px - phi_nx)*twoGD);
        }

        std::array<int, D> order;

        // find largest values of the normal vector
        // (the directions where the normal vector has the smallest deviation from the grid axis)
        if(gradient[0] > gradient[1] ){

            if(gradient[1] > gradient[2]){

                order[0] = 0;
                order[1] = 1;
                order[2] = 2;

            }else{

                if(gradient[0] > gradient[2]){

                    order[0] = 0;
                    order[1] = 2;
                    order[2] = 1;

                }else{

                    order[0] = 2;
                    order[1] = 0;
                    order[2] = 1;
                }
            }
        }else{
            if(gradient[1] <= gradient[2]){

                order[0] = 2;
                order[1] = 1;
                order[2] = 0;
                
            }else{

                if(gradient[2] > gradient[0]){
                    order[0] = 1;
                    order[1] = 2;
                    order[2] = 0;
                }else{
                    order[0] = 1;
                    order[1] = 0;
                    order[2] = 2;
                }

            }
        }

        std::array<T, D> scondOrderDerivatives;

        //calculate all needed derivatives in the planes that have the smalles deviation from the grid axis
        for(int i = 0; i < D; i++){

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = order[i];
            int second_axis;

            if(i == 0){
                second_axis = order[1];
            }else{
                second_axis = order[0];
            }

            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

            //central
            scondOrderDerivatives[i] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*threeGDsq;


        }

        T result = 0.;

        for(int i=0; i<D ; i++)
            result += scondOrderDerivatives[i];

        //mean curvature is trace devided by 2
        return result*0.5;

    }


    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){


        //calculate Gradient

        std::array<T,D> gradient;
        //TODO: think of using upwind normals?
        for(int i = 0; i < D; i++){

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            posUnit[i] = 1;
            negUnit[i] = -1;

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            gradient[i] = std::abs((phi_px - phi_nx)*twoGD);
        }

        std::array<int, D> order;

        // find largest values of the normal vector
        // (the directions where the normal vector has the smallest deviation from the grid axis)
        if(gradient[0] > gradient[1] ){

            if(gradient[1] > gradient[2]){

                order[0] = 0;
                order[1] = 1;
                order[2] = 2;

            }else{

                if(gradient[0] > gradient[2]){

                    order[0] = 0;
                    order[1] = 2;
                    order[2] = 1;

                }else{

                    order[0] = 2;
                    order[1] = 0;
                    order[2] = 1;
                }
            }
        }else{
            if(gradient[1] <= gradient[2]){

                order[0] = 2;
                order[1] = 1;
                order[2] = 0;
                
            }else{

                if(gradient[2] > gradient[0]){
                    order[0] = 1;
                    order[1] = 2;
                    order[2] = 0;
                }else{
                    order[0] = 1;
                    order[1] = 0;
                    order[2] = 2;
                }

            }
        }

        std::array<T, D> scondOrderDerivatives;

        //calculate all needed derivatives in the planes that have the smalles deviation from the grid axis
        for(int i = 0; i < D; i++){

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = order[i];
            int second_axis;

            if(i == 0){
                second_axis = order[1];
            }else{
                second_axis = order[0];
            }

            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

            //central
            scondOrderDerivatives[i] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*threeGDsq;


        }

        T result = 0.;

        for(int i=0; i<D ; i++)
            result += scondOrderDerivatives[i];

        //mean curvature is trace devided by 2
        return result*0.5;

    }




};


template <class T, int D> class variationOfNormals : public baseDerivative<T, D>{

    private:

    //hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator;

    T gridDelta;

    T GD = 0;
    T twoGD = 0;

    public:

    variationOfNormals(T mGD)
    :  gridDelta(mGD){
        GD = 1./gridDelta;
        twoGD = 1./(2.*gridDelta);
    }

    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */

   

    T operator()(hrleCartesianPlaneIterator<hrleDomain<T, D>> & neighborIterator){

        //calculate all needed derivatives in the xy yz and xz plane

        std::array<T, 3> centralDiff;

        std::array<T, 6> oneSidedeDiff;

        std::array<T, 6> derivativesPos;

        std::array<T, 6> derivativesNeg;
        

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_pos  = i;
            int second_pos = (i+1) % D;
            

            posUnit[first_pos] = 1;
            negUnit[first_pos] = -1;
        
            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

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

            posUnit[second_pos] = -1;
            negUnit[second_pos] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();


            //central
            centralDiff[i] = (phi_px - phi_nx)*twoGD;

            //one sided
            oneSidedeDiff[i] = (phi_px - phi_0)*GD;
            oneSidedeDiff[i+3] = (phi_0 - phi_nx)*GD;

            //central outer
            derivativesPos[i]   = (phi_pp - phi_np)*twoGD; 
            derivativesPos[i+3] = (phi_pp - phi_pn)*twoGD; 

            derivativesNeg[i] = (phi_pn- phi_nn)*twoGD; 
            derivativesNeg[i+3] = (phi_np - phi_nn)*twoGD; 

        }

        T n_x = (oneSidedeDiff[0] /
                (std::sqrt((oneSidedeDiff[0]*oneSidedeDiff[0]) + 
                std::pow((derivativesPos[3] + centralDiff[1])*0.5, 2) + 
                std::pow((derivativesPos[5] + centralDiff[2])*0.5, 2)))) 
                -
                (oneSidedeDiff[3] /
                (std::sqrt((oneSidedeDiff[3]*oneSidedeDiff[3]) + 
                std::pow((derivativesNeg[3] + centralDiff[1])*0.5, 2) + 
                std::pow((derivativesNeg[5] + centralDiff[2])*0.5, 2))));

        T n_y = (oneSidedeDiff[1] /
                (std::sqrt((oneSidedeDiff[1]*oneSidedeDiff[1]) + 
                std::pow((derivativesPos[0] + centralDiff[0])*0.5, 2) + 
                std::pow((derivativesPos[1] + centralDiff[2])*0.5, 2)))) 
                -
                (oneSidedeDiff[4] /
                (std::sqrt((oneSidedeDiff[4]*oneSidedeDiff[4]) + 
                std::pow((derivativesNeg[0] + centralDiff[0])*0.5, 2) + 
                std::pow((derivativesNeg[1] + centralDiff[2])*0.5, 2))));

        T n_z = (oneSidedeDiff[2] /
                (std::sqrt((oneSidedeDiff[2]*oneSidedeDiff[2]) + 
                std::pow((derivativesPos[2] + centralDiff[0])*0.5, 2) + 
                std::pow((derivativesPos[4] + centralDiff[1])*0.5, 2)))) 
                -
                (oneSidedeDiff[5] /
                (std::sqrt((oneSidedeDiff[5]*oneSidedeDiff[5]) + 
                std::pow((derivativesNeg[2] + centralDiff[0])*0.5, 2) + 
                std::pow((derivativesNeg[4] + centralDiff[1])*0.5, 2))));

                           
        return (n_x + n_y + n_z);
        

    }

    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){

        //calculate all needed derivatives in the xy yz and xz plane

        std::array<T, 3> centralDiff;

        std::array<T, 6> oneSidedeDiff;

        std::array<T, 6> derivativesPos;

        std::array<T, 6> derivativesNeg;
        

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_pos  = i;
            int second_pos = (i+1) % D;
            

            posUnit[first_pos] = 1;
            negUnit[first_pos] = -1;
        
            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

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

            posUnit[second_pos] = -1;
            negUnit[second_pos] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();


            //central
            centralDiff[i] = (phi_px - phi_nx)*twoGD;

            //one sided
            oneSidedeDiff[i] = (phi_px - phi_0)*GD;
            oneSidedeDiff[i+3] = (phi_0 - phi_nx)*GD;

            //central outer
            derivativesPos[i]   = (phi_pp - phi_np)*twoGD; 
            derivativesPos[i+3] = (phi_pp - phi_pn)*twoGD; 

            derivativesNeg[i] = (phi_pn- phi_nn)*twoGD; 
            derivativesNeg[i+3] = (phi_np - phi_nn)*twoGD; 

        }

        T n_x = (oneSidedeDiff[0] /
                (std::sqrt((oneSidedeDiff[0]*oneSidedeDiff[0]) + 
                std::pow((derivativesPos[3] + centralDiff[1])*0.5, 2) + 
                std::pow((derivativesPos[5] + centralDiff[2])*0.5, 2)))) 
                -
                (oneSidedeDiff[3] /
                (std::sqrt((oneSidedeDiff[3]*oneSidedeDiff[3]) + 
                std::pow((derivativesNeg[3] + centralDiff[1])*0.5, 2) + 
                std::pow((derivativesNeg[5] + centralDiff[2])*0.5, 2))));

        T n_y = (oneSidedeDiff[1] /
                (std::sqrt((oneSidedeDiff[1]*oneSidedeDiff[1]) + 
                std::pow((derivativesPos[0] + centralDiff[0])*0.5, 2) + 
                std::pow((derivativesPos[1] + centralDiff[2])*0.5, 2)))) 
                -
                (oneSidedeDiff[4] /
                (std::sqrt((oneSidedeDiff[4]*oneSidedeDiff[4]) + 
                std::pow((derivativesNeg[0] + centralDiff[0])*0.5, 2) + 
                std::pow((derivativesNeg[1] + centralDiff[2])*0.5, 2))));

        T n_z = (oneSidedeDiff[2] /
                (std::sqrt((oneSidedeDiff[2]*oneSidedeDiff[2]) + 
                std::pow((derivativesPos[2] + centralDiff[0])*0.5, 2) + 
                std::pow((derivativesPos[4] + centralDiff[1])*0.5, 2)))) 
                -
                (oneSidedeDiff[5] /
                (std::sqrt((oneSidedeDiff[5]*oneSidedeDiff[5]) + 
                std::pow((derivativesNeg[2] + centralDiff[0])*0.5, 2) + 
                std::pow((derivativesNeg[4] + centralDiff[1])*0.5, 2))));

                           
        return (n_x + n_y + n_z);
        

    }


};

template <class T, int D> class curvaturGeneralFormulaBigStencil : public baseDerivative<T, D>{

    private:

    T gridDelta;

    T twoGD = 0;
    T fourGDsq = 0;
    T threeGDsq = 0;
    T fourGD = 0;

    public:

    curvaturGeneralFormulaBigStencil(T mGD)
    :  gridDelta(mGD){

        twoGD = 1./(2.*gridDelta);
        fourGD = 1./(4.*gridDelta);
        threeGDsq = 1./(3.*gridDelta*gridDelta);
        fourGDsq = 1./(4.*gridDelta*gridDelta);
    }
    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */


    T operator()(hrleCartesianPlaneIterator<hrleDomain<T, D>> & neighborIterator){

        //array to store derivatives
        std::array<T, 9> d;

        //calculate higher order derivatives

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = i;
            int second_axis = (i+1)%D;

            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();


            //d[i] = (phi_px - phi_nx)*twoGD;

            d[i] = (phi_pp - phi_np + phi_pn - phi_nn)*fourGD;

            d[i+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*threeGDsq;

            d[i+6] = (phi_pp - phi_pn - phi_np + phi_nn)*fourGDsq;

        }

        T norm_grad_pow3 = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;


        //expanded fom of the equation
        return
        //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx + F_yy)
        (d[0]*d[0]*(d[4] + d[5]) + d[1]*d[1]*(d[3] + d[5]) + d[2]*d[2]*(d[3] + d[4]) +

        //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
        -2.*(d[0]*d[1]*d[6] + d[0]*d[2]*d[8] + d[1]*d[2]*d[7]))
                
        // /2*(F_x² + F_y² + F_z²)^(3/2)
        /(2.*norm_grad_pow3);

        //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)

    }

    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){

        //array to store derivatives
        std::array<T, 9> d;

        //calculate higher order derivatives

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = i;
            int second_axis = (i+1)%D;

            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();


            //d[i] = (phi_px - phi_nx)*twoGD;

            d[i] = (phi_pp - phi_np + phi_pn - phi_nn)*fourGD;

            d[i+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*threeGDsq;

            d[i+6] = (phi_pp - phi_pn - phi_np + phi_nn)*fourGDsq;

        }

        T norm_grad_pow3 = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;


        //expanded fom of the equation
        return
        //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx + F_yy)
        (d[0]*d[0]*(d[4] + d[5]) + d[1]*d[1]*(d[3] + d[5]) + d[2]*d[2]*(d[3] + d[4]) +

        //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
        -2.*(d[0]*d[1]*d[6] + d[0]*d[2]*d[8] + d[1]*d[2]*d[7]))
                
        // /2*(F_x² + F_y² + F_z²)^(3/2)
        /(2.*norm_grad_pow3);

        //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)

    }


};

template <class T, int D> class curvaturGeneralFormulaBigStencilBias : public baseDerivative<T, D>{

    private:

    T gridDelta;

    T twoGD = 0;
    T fourGD = 0;
    T fourGDsq = 0;
    T threeGDsq = 0;

    public:

    curvaturGeneralFormulaBigStencilBias(T mGD)
    :  gridDelta(mGD){

        twoGD = 1./(2.*gridDelta);
        fourGD = 1./(4.*gridDelta);
        threeGDsq = 1./(3.*gridDelta*gridDelta);
        fourGDsq = 1./(4.*gridDelta*gridDelta);
    }
    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */

   
    T operator()(hrleCartesianPlaneIterator<hrleDomain<T, D>> & neighborIterator){

        //array to store derivatives
        std::array<T, 9> d;


        //calculate Gradient first for biasing
        std::array<T,D> g;

        for(int i = 0; i < D; i++){

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            posUnit[i] = 1;
            negUnit[i] = -1;

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            //d[i] = (phi_px - phi_nx)*0.5;

            d[i] = (phi_px - phi_nx)*twoGD;

            g[i] = std::abs(d[i]);

        }

                    
        std::array<int, D> order;

        //TODO: konstant die ebenen hinschreiben hinschreiben

        if(g[0] >= g[1] ){

            if(g[1] > g[2]){
                //0,1,2

                order[0] = 0;
                order[1] = 1;
                order[2] = 2;

            }else{

                if(g[0] > g[2]){
                    //0,2,1
                    
                    order[0] = 0;
                    order[1] = 2;
                    order[2] = 1;

                }else{
                    //2,0,1

                    order[0] = 2;
                    order[1] = 0;
                    order[2] = 1;
                }
            }
        }else{
            if(g[1] <= g[2]){
                //2,1,0

                order[0] = 2;
                order[1] = 1;
                order[2] = 0;
                
            }else{

                if(g[2] > g[0]){
                    //1,2,0
                    order[0] = 1;
                    order[1] = 2;
                    order[2] = 0;
                }else{
                    //1,0,2
                    order[0] = 1;
                    order[1] = 0;
                    order[2] = 2;
                }

            }
        }

        //calculate higher order derivatives
       
        //get required ls values

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = order[i]; //i;
            int second_axis;

            if(i == 0){
                second_axis = order[1];
            }else{
                second_axis = order[0];
            }

            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();
           

            // second order derivatives in the same direction
            // order array is needed to put derivatives in the correct position in the derivatives array
            //d[order[i]+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*oneThird;

            d[order[i]+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*threeGDsq;

            // does not improve qulity significantly
            d[order[i]] = (phi_pp - phi_np + phi_pn - phi_nn)*fourGD;

            //d[order[i]+]

        }

        //TODO: check what timig results say
        //For dxdy derivatives the biasing leads to using wrong planes
        //For dxdy derivatives all 3 planes are required
        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = i; //i;
            int second_axis = (i+1)%D;


            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();
          

            d[i+6] = (phi_pp - phi_pn - phi_np + phi_nn)*fourGDsq;

        }


        T norm_grad_pow3 = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;

        //expanded fom of the equation
        return
        //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx + F_yy)
        (d[0]*d[0]*(d[4] + d[5]) + d[1]*d[1]*(d[3] + d[5]) + d[2]*d[2]*(d[3] + d[4]) +

        //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
        -2.*(d[0]*d[1]*d[6] + d[0]*d[2]*d[8] + d[1]*d[2]*d[7]))
                
        // /2*(F_x² + F_y² + F_z²)^(3/2)
        /(2.*norm_grad_pow3);

        //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)

    }

    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){

        //array to store derivatives
        std::array<T, 9> d;


        //calculate Gradient first for biasing
        std::array<T,D> g;

        for(int i = 0; i < D; i++){

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            posUnit[i] = 1;
            negUnit[i] = -1;

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            //d[i] = (phi_px - phi_nx)*0.5;

            d[i] = (phi_px - phi_nx)/(2.*gridDelta);

            g[i] = std::abs(d[i]);

        }

                    
        std::array<int, D> order;

        //TODO: konstant die ebenen hinschreiben hinschreiben

        if(g[0] >= g[1] ){

            if(g[1] > g[2]){
                //0,1,2

                order[0] = 0;
                order[1] = 1;
                order[2] = 2;

            }else{

                if(g[0] > g[2]){
                    //0,2,1
                    
                    order[0] = 0;
                    order[1] = 2;
                    order[2] = 1;

                }else{
                    //2,0,1

                    order[0] = 2;
                    order[1] = 0;
                    order[2] = 1;
                }
            }
        }else{
            if(g[1] <= g[2]){
                //2,1,0

                order[0] = 2;
                order[1] = 1;
                order[2] = 0;
                
            }else{

                if(g[2] > g[0]){
                    //1,2,0
                    order[0] = 1;
                    order[1] = 2;
                    order[2] = 0;
                }else{
                    //1,0,2
                    order[0] = 1;
                    order[1] = 0;
                    order[2] = 2;
                }

            }
        }

        //calculate higher order derivatives
       
        //get required ls values

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = order[i]; //i;
            int second_axis;

            if(i == 0){
                second_axis = order[1];
            }else{
                second_axis = order[0];
            }

            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();
           

            // second order derivatives in the same direction
            // order array is needed to put derivatives in the correct position in the derivatives array
            //d[order[i]+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*oneThird;

            d[order[i]+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*threeGDsq;

            // does not improve qulity significantly
            d[order[i]] = (phi_pp - phi_np + phi_pn - phi_nn)*fourGD;

            //d[order[i]+]

        }

        //TODO: check what timig results say
        //For dxdy derivatives the biasing leads to using wrong planes
        //For dxdy derivatives all 3 planes are required
        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = i; //i;
            int second_axis = (i+1)%D;


            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();
          

            d[i+6] = (phi_pp - phi_pn - phi_np + phi_nn)*fourGDsq;

        }


        T norm_grad_pow3 = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;

        //expanded fom of the equation
        return
        //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx + F_yy)
        (d[0]*d[0]*(d[4] + d[5]) + d[1]*d[1]*(d[3] + d[5]) + d[2]*d[2]*(d[3] + d[4]) +

        //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
        -2.*(d[0]*d[1]*d[6] + d[0]*d[2]*d[8] + d[1]*d[2]*d[7]))
                
        // /2*(F_x² + F_y² + F_z²)^(3/2)
        /(2.*norm_grad_pow3);

        //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)

    }

};

template <class T, int D> class curvaturTest : public baseDerivative<T, D>{

    //TODO: Clean up!

    private:

    T gridDelta;

    T twoGD = 0;
    T fourGD = 0;
    T threeGDsq = 0;
    T fourGDsq = 0;

    public:

    curvaturTest(T mGD)
    :  gridDelta(mGD){
        twoGD = 1./(2.*gridDelta);
        fourGD = 1./(4.*gridDelta);
        threeGDsq = 1./(3.*gridDelta*gridDelta);
        fourGDsq = 1./(4.*gridDelta*gridDelta);
    }

//hrleCartesianPlaneIterator<hrleDomain<T, D>> & neighborIterator

    T operator()(hrleCartesianPlaneIterator<hrleDomain<T, D>> & neighborIterator){

        //array to store derivatives
        std::array<T, 9> d;

        std::array<T, 9> d1;


        //calculate Gradient first for biasing
        std::array<T,D> g;

        for(int i = 0; i < D; i++){

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            posUnit[i] = 1;
            negUnit[i] = -1;

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            d[i] = (phi_px - phi_nx)*twoGD;

            g[i] = std::abs(d[i]);

        }

                    
        std::array<int, D> order;

        //TODO: konstant die ebenen hinschreiben hinschreiben

        if(g[0] > g[1] ){

            if(g[1] > g[2]){
                //0,1,2

                order[0] = 0;
                order[1] = 1;
                order[2] = 2;

            }else{

                if(g[0] > g[2]){
                    //0,2,1
                    
                    order[0] = 0;
                    order[1] = 2;
                    order[2] = 1;

                }else{
                    //2,0,1

                    order[0] = 2;
                    order[1] = 0;
                    order[2] = 1;
                }
            }
        }else{
            if(g[1] <= g[2]){
                //2,1,0

                order[0] = 2;
                order[1] = 1;
                order[2] = 0;
                
            }else{

                if(g[2] > g[0]){
                    //1,2,0
                    order[0] = 1;
                    order[1] = 2;
                    order[2] = 0;
                }else{
                    //1,0,2
                    order[0] = 1;
                    order[1] = 0;
                    order[2] = 2;
                }

            }
        }

        //calculate higher order derivatives
       
        //get required ls values

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = order[i]; //i;
            int second_axis;

            if(i == 0){
                second_axis = order[1];
            }else{
                second_axis = order[0];
            }

            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

            d[order[i]+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*threeGDsq;

            d[order[i]] = (phi_pp - phi_np + phi_pn - phi_nn)*fourGD;

        }

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = order[i]; //i;
            int second_axis;

            if(i == 2){
                second_axis = order[1];
            }else{
                second_axis = order[2];
            }

            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

            d1[order[i]+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*threeGDsq;

            d1[order[i]] = (phi_pp - phi_np + phi_pn - phi_nn)*fourGD;

        }

        //TODO: check what timig results say
        //For dxdy derivatives the biasing leads to using wrong planes
        //For dxdy derivatives all 3 planes are required
        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = i; //i;
            int second_axis = (i+1)%D;


            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();
          

            d[i+6] = (phi_pp - phi_pn - phi_np + phi_nn)*fourGDsq;

            d1[i+6] = (phi_pp - phi_pn - phi_np + phi_nn)*fourGDsq;

        }


        T norm_grad_pow3 = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;

                //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx + F_yy)
        T erg1 = (d[0]*d[0]*(d[4] + d[5]) + d[1]*d[1]*(d[3] + d[5]) + d[2]*d[2]*(d[3] + d[4]) +

        //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
        -2.*(d[0]*d[1]*d[6] + d[0]*d[2]*d[8] + d[1]*d[2]*d[7]))
                
        // /2*(F_x² + F_y² + F_z²)^(3/2)
        /(2.*norm_grad_pow3);

        norm_grad_pow3 = std::sqrt(d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2]);
        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;

                //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx + F_yy)
        T erg2 = (d1[0]*d1[0]*(d1[4] + d1[5]) + d1[1]*d1[1]*(d1[3] + d1[5]) + d1[2]*d1[2]*(d1[3] + d1[4]) +

        //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
        -2.*(d1[0]*d1[1]*d1[6] + d1[0]*d1[2]*d1[8] + d1[1]*d1[2]*d1[7]))
                
        // /2*(F_x² + F_y² + F_z²)^(3/2)
        /(2.*norm_grad_pow3);

        //expanded fom of the equation
        return (erg1 +erg2)/2.;



        //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)

    }

    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){

        //array to store derivatives
        std::array<T, 9> d;

        std::array<T, 9> d1;


        //calculate Gradient first for biasing
        std::array<T,D> g;

        for(int i = 0; i < D; i++){

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            posUnit[i] = 1;
            negUnit[i] = -1;

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            d[i] = (phi_px - phi_nx)*twoGD;

            g[i] = std::abs(d[i]);

        }

                    
        std::array<int, D> order;

        //TODO: konstant die ebenen hinschreiben hinschreiben

        if(g[0] > g[1] ){

            if(g[1] > g[2]){
                //0,1,2

                order[0] = 0;
                order[1] = 1;
                order[2] = 2;

            }else{

                if(g[0] > g[2]){
                    //0,2,1
                    
                    order[0] = 0;
                    order[1] = 2;
                    order[2] = 1;

                }else{
                    //2,0,1

                    order[0] = 2;
                    order[1] = 0;
                    order[2] = 1;
                }
            }
        }else{
            if(g[1] <= g[2]){
                //2,1,0

                order[0] = 2;
                order[1] = 1;
                order[2] = 0;
                
            }else{

                if(g[2] > g[0]){
                    //1,2,0
                    order[0] = 1;
                    order[1] = 2;
                    order[2] = 0;
                }else{
                    //1,0,2
                    order[0] = 1;
                    order[1] = 0;
                    order[2] = 2;
                }

            }
        }

        //calculate higher order derivatives
       
        //get required ls values

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = order[i]; //i;
            int second_axis;

            if(i == 0){
                second_axis = order[1];
            }else{
                second_axis = order[0];
            }

            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

            d[order[i]+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*threeGDsq;

            d[order[i]] = (phi_pp - phi_np + phi_pn - phi_nn)*fourGD;

        }

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = order[i]; //i;
            int second_axis;

            if(i == 2){
                second_axis = order[1];
            }else{
                second_axis = order[2];
            }

            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

            d1[order[i]+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*threeGDsq;

            d1[order[i]] = (phi_pp - phi_np + phi_pn - phi_nn)*fourGD;

        }

        //TODO: check what timig results say
        //For dxdy derivatives the biasing leads to using wrong planes
        //For dxdy derivatives all 3 planes are required
        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = i; //i;
            int second_axis = (i+1)%D;


            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = 1;
            negUnit[second_axis] = 1;

            T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_axis] = -1;
            negUnit[second_axis] = -1;

            T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[first_axis] = 0;
            negUnit[first_axis] = 0;

            posUnit[second_axis] = 1;
            negUnit[second_axis] = -1;

            T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();
          

            d[i+6] = (phi_pp - phi_pn - phi_np + phi_nn)*fourGDsq;

            d1[i+6] = (phi_pp - phi_pn - phi_np + phi_nn)*fourGDsq;

        }


        T norm_grad_pow3 = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;

                //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx + F_yy)
        T erg1 = (d[0]*d[0]*(d[4] + d[5]) + d[1]*d[1]*(d[3] + d[5]) + d[2]*d[2]*(d[3] + d[4]) +

        //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
        -2.*(d[0]*d[1]*d[6] + d[0]*d[2]*d[8] + d[1]*d[2]*d[7]))
                
        // /2*(F_x² + F_y² + F_z²)^(3/2)
        /(2.*norm_grad_pow3);

        norm_grad_pow3 = std::sqrt(d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2]);
        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;

                //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx + F_yy)
        T erg2 = (d1[0]*d1[0]*(d1[4] + d1[5]) + d1[1]*d1[1]*(d1[3] + d1[5]) + d1[2]*d1[2]*(d1[3] + d1[4]) +

        //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
        -2.*(d1[0]*d1[1]*d1[6] + d1[0]*d1[2]*d1[8] + d1[1]*d1[2]*d1[7]))
                
        // /2*(F_x² + F_y² + F_z²)^(3/2)
        /(2.*norm_grad_pow3);

        //expanded fom of the equation
        return (erg1 +erg2)/2.;



        //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)

    }




};

