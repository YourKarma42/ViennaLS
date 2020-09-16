#include <iostream>

#include <lsDomain.hpp>

#include <hrleSparseBoxIterator.hpp>

#include "hrleTestIterator.hpp"

//TODO: calculate grid delta values in consturctor



template <class T, int D> class curvaturGeneralFormula{

    private:

    //hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator;

    T gridDelta;

    public:

    curvaturGeneralFormula(T mGD)
    :  gridDelta(mGD){
    }
    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */
    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){
    //T operator()(hrleCartesianPlaneIterator<hrleDomain<T, D>> & neighborIterator){

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

            //T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            //T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();
           
/*
            // first order derivative
            d[i] = (phi_px - phi_nx)*0.5;

            // second order derivatives in the same direction
            d[i+3] = (phi_px - 2.*phi_0 + phi_nx);

            d[i+6] = (phi_pp - phi_pn -phi_np + phi_nn)*0.25;
*/


            // first order derivative
            d[i] = (phi_px - phi_nx)/(2.*gridDelta);

            // second order derivatives in the same direction
            d[i+3] = (phi_px - 2.*phi_0 + phi_nx)/(gridDelta*gridDelta);

            d[i+6] = (phi_pp - phi_pn -phi_np + phi_nn)/(4.*gridDelta*gridDelta);

        }


        T norm_grad_pow3 = 0.;
        for(int i = 0; i < D; i++){

            norm_grad_pow3 += d[i]*d[i];
        }
        norm_grad_pow3 = std::sqrt(norm_grad_pow3);

        //d[0] = d[0]/norm_grad_pow3;        
        //d[1] = d[1]/norm_grad_pow3;

        //norm_grad_pow3 = 1.;

        //std::cout << norm_grad_pow3 << std::endl;

        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;



        //TODO: not a clean solution think of something different
/*
        if((d[3]*d[1]*d[1] - 2.*d[1]*d[0]*d[6] + d[4]*d[0]*d[0])/(norm_grad_pow3) < 0.024){
            std::cout << d[3] << std::endl;
            std::cout << d[4] << std::endl;
            std::cout << norm_grad_pow3 << std::endl;
            std::cout << "mist" << std::endl;
        }
        */

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
    T operator()(hrleCartesianPlaneIterator<hrleDomain<T, D>> & neighborIterator){

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

            //T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
            //T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();
           
/*
            // first order derivative
            d[i] = (phi_px - phi_nx)*0.5;

            // second order derivatives in the same direction
            d[i+3] = (phi_px - 2.*phi_0 + phi_nx);

            d[i+6] = (phi_pp - phi_pn -phi_np + phi_nn)*0.25;
*/


            // first order derivative
            d[i] = (phi_px - phi_nx)/(2.*gridDelta);

            // second order derivatives in the same direction
            d[i+3] = (phi_px - 2.*phi_0 + phi_nx)/(gridDelta*gridDelta);

            d[i+6] = (phi_pp - phi_pn -phi_np + phi_nn)/(4.*gridDelta*gridDelta);

        }


        T norm_grad_pow3 = 0.;
        for(int i = 0; i < D; i++){

            norm_grad_pow3 += d[i]*d[i];
        }
        norm_grad_pow3 = std::sqrt(norm_grad_pow3);

        //d[0] = d[0]/norm_grad_pow3;        
        //d[1] = d[1]/norm_grad_pow3;

        //norm_grad_pow3 = 1.;

        //std::cout << norm_grad_pow3 << std::endl;

        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;



        //TODO: not a clean solution think of something different
/*
        if((d[3]*d[1]*d[1] - 2.*d[1]*d[0]*d[6] + d[4]*d[0]*d[0])/(norm_grad_pow3) < 0.024){
            std::cout << d[3] << std::endl;
            std::cout << d[4] << std::endl;
            std::cout << norm_grad_pow3 << std::endl;
            std::cout << "mist" << std::endl;
        }
        */

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


template <class T, int D> class curvaturShapeDerivatives1{

    private:

    //hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator;

    T gridDelta;

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

    T operator()(hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType> & neighborIterator){

        std::array<T, 3> scondOrderDerivatives;
       
        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            //hrleVectorType<hrleIndexType, D> test = neighborIterator.getIndices();

            //TODO: ask people if there is a more elegant solution
            //int first_pos  = i;
            //int second_pos = (i+1) % D;
        
            //get required ls values
            T phi_0 = neighborIterator.getCenter().getValue();

            T phi_px = neighborIterator.getNeighbor(i).getValue();
            T phi_nx = neighborIterator.getNeighbor(i+D).getValue();


            //central
            //scondOrderDerivatives[i] = (phi_px - 2.*phi_0 + phi_nx);  

            scondOrderDerivatives[i] = (phi_px - 2.*phi_0 + phi_nx)/(gridDelta*gridDelta);  
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

            //hrleVectorType<hrleIndexType, D> test = neighborIterator.getIndices();

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
            scondOrderDerivatives[i] = (phi_px - 2.*phi_0 + phi_nx);  

        }


        T result = 0.;

        for(int i=0; i<D ; i++)
            result += scondOrderDerivatives[i];
        //mean curvature is trace devided by 2
        return result*0.5;

    }

};



template <class T, int D> class curvaturShapeDerivatives2{

    private:

    T gridDelta;

    const T oneThird = 1./3.;

    public:

    curvaturShapeDerivatives2(T mGD)
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

        std::array<T, D> scondOrderDerivatives;

        for(int i = 0; i < D; i ++){
       

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            //hrleVectorType<hrleIndexType, D> test = neighborIterator.getIndices();

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
            scondOrderDerivatives[i] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)/(3*gridDelta*gridDelta);
        }

        T result = 0.;

        for(int i=0; i<D ; i++)
            result += scondOrderDerivatives[i];
        //mean curvature is trace devided by 2
        return result*0.5;

    }

};

template <class T, int D> class curvaturShapeBias{

    private:

    T gridDelta;

    const T oneThird = 1./3.;

    public:

    curvaturShapeBias(T mGD)
    :  gridDelta(mGD){
    }
    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */

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

            gradient[i] = std::abs((phi_px - phi_nx)/(2*gridDelta));
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

            //TODO: Explain this part if it works!
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
            scondOrderDerivatives[i] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)/(3.*gridDelta*gridDelta);


        }

        T result = 0.;

        for(int i=0; i<D ; i++)
            result += scondOrderDerivatives[i];

        //mean curvature is trace devided by 2
        return result*0.5;

    }




};


template <class T, int D> class variationOfNormals{

    private:

    //hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator;

    T gridDelta;

    const T oneThird = 1./3.;

    public:

    variationOfNormals(T mGD)
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

            //hrleVectorType<hrleIndexType, D> test = neighborIterator.getIndices();

            //TODO: ask people if there is a more elegant solution
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
            centralDiff[i] = (phi_px - phi_nx)*0.5;

            //one sided
            oneSidedeDiff[i] = (phi_px - phi_0);
            oneSidedeDiff[i+3] = (phi_0 - phi_nx);

            //central outer
            derivativesPos[i]   = (phi_pp - phi_np)*0.5; 
            derivativesPos[i+3] = (phi_pp - phi_pn)*0.5; 

            derivativesNeg[i] = (phi_pn- phi_nn)*0.5; 
            derivativesNeg[i+3] = (phi_np - phi_nn)*0.5; 

        }

        //test if direct multiplication is faster   
            //n_x_pos - n_x_neg
        T n_x = (oneSidedeDiff[0] /
                (2*std::sqrt((oneSidedeDiff[0]*oneSidedeDiff[0]) + 
                std::pow((derivativesPos[3] + centralDiff[1])*0.5, 2) + 
                std::pow((derivativesPos[5] + centralDiff[2])*0.5, 2)))) 
                -
                (oneSidedeDiff[3] /
                (2*std::sqrt((oneSidedeDiff[3]*oneSidedeDiff[3]) + 
                std::pow((derivativesNeg[3] + centralDiff[1])*0.5, 2) + 
                std::pow((derivativesNeg[5] + centralDiff[2])*0.5, 2))));

        T n_y = (oneSidedeDiff[1] /
                (2*std::sqrt((oneSidedeDiff[1]*oneSidedeDiff[1]) + 
                std::pow((derivativesPos[0] + centralDiff[0])*0.5, 2) + 
                std::pow((derivativesPos[1] + centralDiff[2])*0.5, 2)))) 
                -
                (oneSidedeDiff[4] /
                (2*std::sqrt((oneSidedeDiff[4]*oneSidedeDiff[4]) + 
                std::pow((derivativesNeg[0] + centralDiff[0])*0.5, 2) + 
                std::pow((derivativesNeg[1] + centralDiff[2])*0.5, 2))));

        T n_z = (oneSidedeDiff[2] /
                (2*std::sqrt((oneSidedeDiff[2]*oneSidedeDiff[2]) + 
                std::pow((derivativesPos[2] + centralDiff[0])*0.5, 2) + 
                std::pow((derivativesPos[4] + centralDiff[1])*0.5, 2)))) 
                -
                (oneSidedeDiff[5] /
                (2*std::sqrt((oneSidedeDiff[5]*oneSidedeDiff[5]) + 
                std::pow((derivativesNeg[2] + centralDiff[0])*0.5, 2) + 
                std::pow((derivativesNeg[4] + centralDiff[1])*0.5, 2))));

                           
        return n_x + n_y + n_z;
        

    }


};

template <class T, int D> class curvaturGeneralFormulaBigStencil{

    private:

    T gridDelta;

    const T oneThird = 1./3.;

    public:

    curvaturGeneralFormulaBigStencil(T mGD)
    :  gridDelta(mGD){
    }
    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */

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
           

            // second order derivatives in the same direction
            // order array is needed to put derivatives in the correct position in the derivatives array

            //d[i] = (phi_px - phi_nx)*0.5;

            //d[i+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*oneThird;

            //d[i+6] = (phi_pp - phi_pn - phi_np + phi_nn)*0.25;


            d[i] = (phi_px - phi_nx)/(2.*gridDelta);

            d[i+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)/(3.*gridDelta*gridDelta);

            d[i+6] = (phi_pp - phi_pn - phi_np + phi_nn)/(4.*gridDelta*gridDelta);

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

template <class T, int D> class curvaturGeneralFormulaBigStencilBias{

    private:

    T gridDelta;

    const T oneThird = 1./3.;

    public:

    curvaturGeneralFormulaBigStencilBias(T mGD)
    :  gridDelta(mGD){
    }
    /*

    Slice of a stencil: x axis horizontal y axis vertical

        phi_np | phi_py | phi_pp
        phi_nx | phi_0  | phi px
        phi_nn | phi_ny | phi_nn
    */

    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){

        //array to store derivatives
        std::array<T, 9> d;


        //calculate Gradient first for biasing
        std::array<T,D> g;

        //TODO: think of using upwind normals?
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
           

            // second order derivatives in the same direction
            // order array is needed to put derivatives in the correct position in the derivatives array
            //d[order[i]+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)*oneThird;

            d[order[i]+3] = (phi_pp - 2.*phi_py + phi_np + phi_px -2.*phi_0 + phi_nx + phi_pn - 2.*phi_ny + phi_nn)/(3.*gridDelta*gridDelta);

            // does not improve qulity significantly
            d[order[i]] = (phi_pp - phi_np + phi_pn - phi_nn)/(4.*gridDelta);

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
           

            d[i+6] = (phi_pp - phi_pn - phi_np + phi_nn)/(4.*gridDelta*gridDelta);



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

template <class T, int D> class curvaturTest{

    //TODO: Clean up!

    private:

    T gridDelta;

    const T oneThird = 1./3.;

    const T threeThirteen = 3./13.;

    public:

    curvaturTest(T mGD)
    :  gridDelta(mGD){
        //TODO:check width of levelset! print error if not big enough
    }

    //TODO:reference to paper
    //only for 3D
    T operator()(hrleSparseBoxIterator<hrleDomain<T, D>> & neighborIterator){

        //Center of the stencil is needed in each sum
        //the center also only occures times 2
        T twoPhi_0 = 2.*neighborIterator.getCenter().getValue();

        //calculate corner derivatives (diagonals) (the same for each direction)

        T cornerDerivatives = 0.;

        cornerDerivatives +=  neighborIterator.getNeighbor(hrleVectorType<hrleIndexType, D>(1,1,1)).getValue() + 
                                neighborIterator.getNeighbor(hrleVectorType<hrleIndexType, D>(-1,-1,-1)).getValue() - twoPhi_0;

        cornerDerivatives +=  neighborIterator.getNeighbor(hrleVectorType<hrleIndexType, D> (-1,1,1)).getValue() + 
                                neighborIterator.getNeighbor(hrleVectorType<hrleIndexType, D> (1,-1,-1)).getValue() - twoPhi_0;

        cornerDerivatives +=  neighborIterator.getNeighbor(hrleVectorType<hrleIndexType, D>(-1,-1,1)).getValue() + 
                                neighborIterator.getNeighbor(hrleVectorType<hrleIndexType, D>(1,1,-1)).getValue() - twoPhi_0;

        cornerDerivatives +=  neighborIterator.getNeighbor(hrleVectorType<hrleIndexType, D>(1,-1,1)).getValue() + 
                                neighborIterator.getNeighbor(hrleVectorType<hrleIndexType, D>(-1,1,-1)).getValue() - twoPhi_0;                        

        cornerDerivatives = cornerDerivatives;

        //Calculate the derivates in the plane and in center of the stencil
        std::array<T, 2*D> planeDerivatives;

        std::array<T, D> centerDerivatives;

        for (int i = 0; i < D; i++) {

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            int first_axis = i;
            int second_axis = (i+1)%D;

            posUnit[first_axis] = 1;
            negUnit[first_axis] = -1;

            //get required ls values
            //T phi_0 = neighborIterator.getCenter().getValue();

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

            centerDerivatives[i] = (phi_px - twoPhi_0 + phi_nx);  

            planeDerivatives[i] = (phi_pp - twoPhi_0 + phi_nn);  

            planeDerivatives[i+D] = (phi_np - twoPhi_0 + phi_pn);  

        }
        
        //calculate the actual second order derivatives
        std::array<T, D> scondOrderDerivatives;



        scondOrderDerivatives[0] =  threeThirteen * (centerDerivatives[0] + 
                                    (planeDerivatives[0] + planeDerivatives[0+D]) * 0.5 
                                    + cornerDerivatives * oneThird);

        scondOrderDerivatives[1] =  threeThirteen * (centerDerivatives[1] + 
                                    (planeDerivatives[1] + planeDerivatives[1+D]) * 0.5 
                                    + cornerDerivatives * oneThird);   

        scondOrderDerivatives[2] =  threeThirteen * (centerDerivatives[2] + 
                                    (planeDerivatives[2] + planeDerivatives[2+D]) * 0.5 
                                    + cornerDerivatives * oneThird);

        return (threeThirteen*((centerDerivatives[0] + centerDerivatives[1] + centerDerivatives[2])
                + 0.5*(planeDerivatives[0] + planeDerivatives[0+D] + planeDerivatives[1] + planeDerivatives[1+D] + planeDerivatives[2] + planeDerivatives[2+D])
                + oneThird*(cornerDerivatives) ))
        
        
        *0.5;// *threeThirteen *+ 
                //0.5*(planeDerivatives[0]+planeDerivatives[1]+planeDerivatives[2] +planeDerivatives[0+D]+planeDerivatives[1+D]+planeDerivatives[2+D]) + 
                //cornerDerivatives * oneThird;


        T result = 0.;

        for(int i=0; i<D ; i++)
            result += scondOrderDerivatives[i];
        //mean curvature is trace devided by 2
        return result*0.5;

    }




};

