#pragma once
#include <memory>
#include "discretization/discretization.h"
#include <cmath>

class PressureSolver
{
    public:
        
        PressureSolver (std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

        virtual void solve();

    protected:
        void setBoundaryValues();
        void computeResiduum();
        std::shared_ptr< Discretization > 	discretization_           ;
        double 	                            epsilon_                  ;
        int 	                            maximumNumberOfIterations_;
        double                              residuum_                 ;
        int                           N_;
        };
