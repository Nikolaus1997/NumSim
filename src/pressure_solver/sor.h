#pragma once
#include "pressure_solver/pressure_solver.h"

/**
 * This class overrides the solve method from Discretization with the SOR method.
*/

class SOR : public PressureSolver
{
    public:
        //constructor
        SOR (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations, double omega);
    
        //solve the system of the Poisson equation for pressure using SOR
        void solve();
    
    private:
        double omega_;
};
