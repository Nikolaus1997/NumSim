#pragma once
#include <memory>
#include "discretization/discretization.h"

class PressureSolver
{
    public:
        
        PressureSolver (std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

        virtual void solve();

    protected:
        void setBoundaryValues();
        
        std::shared_ptr< Discretization > 	discretization_           ;
        double 	                            epsilon_                  ;
        int 	                            maximumNumberOfIterations_;
};
