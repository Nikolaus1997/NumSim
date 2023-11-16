#pragma once
#include <memory>
#include "discretization/discretization.h"
#include <cmath>

/**
 * This class provides functions and other parameters for the solution of the Poisson problem for the pressure.
*/

class PressureSolver
{
    public:
        //constructor
        PressureSolver (std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

        //solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid 
        virtual void solve();

    protected:
        //set boundary values as homogeneous Neumann boundary values
        void setBoundaryValues();
        void computeResiduum();
        std::shared_ptr< Discretization > 	discretization_           ;
        double 	                            epsilon_                  ;
        int 	                            maximumNumberOfIterations_;
        double                              residuum_                 ;
        int                           N_;
        };
