#ifndef LS_CALCULATE_CURVATURES_HPP
#define LS_CALCULATE_CURVATURES_HPP

#include <lsPreCompileMacros.hpp>

#include <algorithm>

#include <hrleSparseBoxIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain.hpp>
#include <lsMessage.hpp>


//TODO: nen namespace machen?

//TODO describe class!!!!

//TODO: At the moment only one lvl set function

//TODO: Neighbour iterator komplett überarbeiten

//TODO:create max value for level set value

template <class T, int D> class lsCalculateCurvatures {

private:

    bool first = false;

    lsDomain<T, D> &levelSet;

    //TODO::different solution for this iterator
    //Inf values occur when the level set is not expanded previously
    hrleSparseBoxIterator<hrleDomain<T, D>> neighborIterator;

    T denTwoDelta, denDeltaSquared, denFourDeltaSquared; 

    std::vector<T> gaussianCurvatures;

    std::vector<T> meanCurvatures;

    //TODO:remove just for testing normal vectors

    std::vector<std::array<T, D>> normals;

               

    //Calculates all partial derivates of order 2 and returns them in a vector of the form
    //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)
    //All partial derivatives of order 2 are equivalent due to schwarzes rule (F_xy = F_yx) 
    std::vector<T> calc_partial_derivatives(const hrleVectorType<hrleIndexType, D> &indices){
    //int calc_partial_derivatives(const hrleVectorType<hrleIndexType, D> &indices){


        //maby private and only overwrite?
        std::vector<T> derivatives(9);

        double gridDelta = levelSet.getGrid().getGridDelta();

        // move neighborIterator to current position
        neighborIterator.goToIndicesSequential(indices);

        //std::cout << "vals from loop:"  << std::endl;

        //This calculation strategy only works for 3 dimensions!
        for (int i = 0; i < D; i++) { 

            hrleVectorType<hrleIndexType, D> posUnit(0);
            hrleVectorType<hrleIndexType, D> negUnit(0);

            posUnit[i] = 1;
            negUnit[i] = -1;

            int second_pos = (i+1) % D;

            //get required ls values
            const T phi_0 = neighborIterator.getCenter().getValue();

            const T phi_p = neighborIterator.getNeighbor(posUnit).getValue();
            const T phi_n = neighborIterator.getNeighbor(negUnit).getValue();
            posUnit[second_pos] = 1;
            negUnit[second_pos] = 1;

            const T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
            const T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

            posUnit[second_pos] = -1;
            negUnit[second_pos] = -1;

            const T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
            const T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

            //This implementation currently uses the central difference quotient
            //TODO: TRY to avoid division
            //TODO: see if we need division by grid delta (multiplication is faster!)
            //TODO: deviding by griddelta dosnt change during the whole algorithm so calc it once than multiply
            //TODO: Grid delta only needed when using real phi values
            //Calculate partial derivatives for gradient and hessian
            //Calc F_i

            //derivatives[i] = ((phi_p - phi_0 ) + (phi_0 - phi_n))*0.5;

            derivatives[i] = (phi_p - phi_n)*0.5;   // /(2*gridDelta);

            //derivatives[i] = (phi_p - phi_n)*denTwoDelta;

            //Calc F_ii

            derivatives[i+3] = (phi_p - 2*phi_0 + phi_n); // /(gridDelta*gridDelta);


            //derivatives[i+3] = (phi_p - 2*phi_0 + phi_n)*denDeltaSquared;

            //Calc F_i(i mod D-1)
            
            derivatives[i+6] = (phi_pp - phi_pn -phi_np + phi_nn)*0.5; // /(2*gridDelta*gridDelta);

            //T test = (phi_pp - phi_p - phi_);

            //test = test;

            //derivatives[i+6] = (phi_pp - phi_pn -phi_np + phi_nn)*denFourDeltaSquared;<

            //alternative formulaiton
            //See wiki need x-y plane in one!

            //stuff für debug_________________________________________________________________________________________
           
            /*std::cout << "phi 0    phi n    phi p" << std::endl;
            std::cout << phi_0 << " " << phi_n << " " << phi_p << " " << std::endl;
            std::cout << "phi nn    phi pp" << std::endl;
            std::cout << phi_nn << " " << phi_pp << " " << std::endl;
            std::cout << "phi np    phi pn" << std::endl;
            std::cout << phi_np << " " << phi_pp << " " << std::endl;

            int tmp;
            std::cin >> tmp;*/
            
            //stuff für debug_________________________________________________________________________________________

        }

       

        //stuff für debug_________________________________________________________________________________________
       /* int count = 0;

         std::cout <<  std::endl;

        std::cout << "normal iterator" << std::endl;

        for(int i = 0; i < 27; i++){

           // if (neighborIterator.getNeighbor(i).getValue() > 5.0){
           //     count++;
                std::cout << "i: " << i << "  " << neighborIterator.getNeighbor(i).getValue() << std::endl; 

                //;neighborIterator.getNeighbor(i).getValue()
           // }
        }

        //for(auto v : derivatives)
        //    std::cout << v << " ";
            
        std::cout << std::endl;

        std::cin >> count;*/

        //stuff für debug_________________________________________________________________________________________

        return derivatives;
       

    }


    T calc_gaussian_curvature(std::vector<T> d){

        //TODO: write used formula and paper with reference


        T norm_grad_sqared = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
        norm_grad_sqared = norm_grad_sqared*norm_grad_sqared;

        //expanded fom of the equation
        return
        //  - (F_x²(F_yz²-F_yyF_zz)   +     F_y²(F_zx²-F_xxF_zz)     +     F_z²(F_xy²-F_xxFyy) +
        -(d[0]*d[0]*(d[7]*d[7]-d[4]*d[5]) + d[1]*d[1]*(d[8]*d[8]-d[3]*d[5]) + d[2]*d[2]*(d[6]*d[6]-d[3]*d[4]) +

        //2*[F_xF_y(2F_zzF_xy-2F_zxF_yz) +
        2*(d[0]*d[1]*(d[5]*d[6]-d[8]*d[7]) +

        //F_xF_z(2F_yyF_zx-2F_xyF_yz) +
        d[0]*d[2]*(d[4]*d[8]-d[6]*d[7]) +

        //F_yF_z(2F_xxF_yz-2F_xyF_zx))
        d[1]*d[2]*(d[3]*d[7]-d[6]*d[8])))
        
        // /(F_x² + F_y² + F_z²)²
        /norm_grad_sqared;

         //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)
    }

    T calc_mean_curvature(std::vector<T> d){


        //TODO: write used formula and paper with reference


        T norm_grad_pow3 = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        norm_grad_pow3 = norm_grad_pow3*norm_grad_pow3*norm_grad_pow3;

        //expanded fom of the equation
        return
        //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx - F_yy)
        (d[0]*d[0]*(d[4] + d[5]) + d[1]*d[1]*(d[3] + d[5]) + d[2]*d[2]*(d[3] + d[4]) 

        //-2*[F_xF_yF_xy +
        -2*(d[0]*d[1]*d[6] +

        //F_xF_zF_xz +
        d[0]*d[2]*d[8] +

        //F_yF_zF_yz
        d[1]*d[2]*d[7]))
        
        // /(F_x² + F_y² + F_z²)^(3/2)
        /norm_grad_pow3;

         //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)

    }

    //TODO:remove just for testing normal vectors
    std::array<T, D> tmp_calc_normals(std::vector<T> d){

        std::array<T, D> n;

        T norm_grad = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

        n[0] = d[0]/norm_grad;

        n[1] = d[1]/norm_grad;

        n[2] = d[2]/norm_grad;

        /*for(auto v: n){
            std::cout << v << " ";
        }
        std::cout << std::endl;

        int tmp;
        std::cin >> tmp;*/
        

        return n;

    }





public:

     lsCalculateCurvatures(lsDomain<T, D> &passedlsDomain)
    : levelSet(passedlsDomain), neighborIterator(levelSet.getDomain(), 1) {
        
        //for second dreivatives a bigger stencil is needed so expand the level-set
        //TODO: expand levelset here and create iterator!
        
        //precalculate denominaters of fractions for finite differences
        //TODO: not necessary maby compiler optimizes it away?
        double gridDelta = levelSet.getGrid().getGridDelta();

        denTwoDelta = 1/(2*gridDelta);
        
        denDeltaSquared = 1/(gridDelta*gridDelta);  
        
        denFourDeltaSquared = 1/(4*gridDelta*gridDelta);

    }

    //TODO: remove dont need
    static void prepareLS(lsDomain<T, D> &passedlsDomain) {
        int order = 1;
        lsExpand<T, D>(passedlsDomain, 2 * (order + 2) + 1).apply();
    }

    std::vector<T>& getGaussianCurvature(){ 

        //TODO: implement for 2 dim case und dann die andere curvature zurückgeben

        if(D < 3){
            lsMessage::getInstance()
            .addWarning("Only implicit fucntions (Surfaces) of dimension >= 3 have a Gaussian Curvature! Please use the approptiate function!")
            .print();
            return meanCurvatures;
        }
        
        return gaussianCurvatures; 
    }

    std::vector<T>& getMeanCurvature(){ 

        //TODO: implement for 2 dim case und dann die andere curvature zurückgeben

        if(D < 3){
            lsMessage::getInstance()
            .addWarning("Only implicit fucntions (Surfaces) of dimension >= 3 have a Mean curvature! Please use the approptiate function!")
            .print();
            return meanCurvatures;
        }
        
        return meanCurvatures; 
    }

    std::vector<std::array<T, D>>& getNormals(){ 

        //TODO: implement for 2 dim case und dann die andere curvature zurückgeben

        if(D < 3){
            lsMessage::getInstance()
            .addWarning("Only implicit fucntions (Surfaces) of dimension >= 3 have a Mean curvature! Please use the approptiate function!")
            .print();
            //return meanCurvatures;
        }
        
        return normals; 
    }



    



    
    void apply(){

        //TODO: Try to expand levelset here!

        auto &topDomain = levelSet.getDomain();
        auto &grid = levelSet.getGrid();

        //TODO: DONT UNDERSTAND
        double pointsPerSegment =
        double(2 * levelSet.getDomain().getNumberOfPoints()) /
        double(levelSet.getLevelSetWidth());

        //TODO: maby do in consturctor?
        gaussianCurvatures.reserve(pointsPerSegment);
        meanCurvatures.reserve(pointsPerSegment);
        normals.reserve(pointsPerSegment);

        int p = 0;

        //DONT UNDERSTAND
        hrleVectorType<hrleIndexType, D> startVector =
        (p == 0) ? grid.getMinGridPoint()
             : topDomain.getSegmentation()[p - 1];

        //DONT UNDERSTAND
        hrleVectorType<hrleIndexType, D> endVector =
        (p != static_cast<int>(topDomain.getNumberOfSegments() - 1))
            ? topDomain.getSegmentation()[p]
            : grid.incrementIndices(grid.getMaxGridPoint());

        if(D == 3){


            for (hrleSparseIterator<typename lsDomain<T, D>::DomainType> it(
                topDomain, startVector);
                it.getStartIndices() < endVector; ++it) {


                if (!it.isDefined()) {
                    continue;
                } else if (std::abs(it.getValue()) > 0.5) {
                    // push an empty vector to keep ordering correct
                    std::array<T, D> tmp = {};
                    normals.push_back(tmp);
                    //TODO: think of good return value maby inf?
                    gaussianCurvatures.push_back(1337.1337);
                    meanCurvatures.push_back(1337.1337);
                    continue;
                }

                //TODO: Remove
                //checks if the run is defined
                /*if (!it.isDefined() || std::abs(it.getValue()) > 0.5){
                    continue;
                }else{

                }*/

                //make method member
                std::vector<T> derivatives = calc_partial_derivatives(it.getStartIndices());

                gaussianCurvatures.push_back(calc_gaussian_curvature(derivatives));

                meanCurvatures.push_back(calc_mean_curvature(derivatives));

                normals.push_back(tmp_calc_normals(derivatives));

                //TODO: REMOVE
                if(std::isnan(gaussianCurvatures.back()))
                    std::cout << "error!!!!!!!!! ";

                // double test = calc_gaussian_curvature(derivatives);

                //std::cout << "test " << test << std::endl;

               /* std::cout << "derivatives: " << std::endl;

                for(auto val: derivatives){
                    std::cout << (double)val << " " << std::endl;

                }

                std::cout << std::endl;

                std::cout << "Gauss: " << gaussianCurvatures.back() << std::endl;

                

                std::cin >> p;
                */

            }

        }else if(D == 2){
            lsMessage::getInstance()
                .addWarning("Curvature calculation for 2D not implemented yet!")
                .print();

        }else{
            lsMessage::getInstance()
                .addWarning("Curvature Calculation is only implemented for D = 2 and D = 3.")
                .print();
            
        }


    }


};



// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsCalculateCurvatures)

#endif // LS_CALCULATE_CURVATURES_HPP