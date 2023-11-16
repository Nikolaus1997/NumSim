#pragma once
#include "pressure_solver/pressure_solver.h"
#include <cmath>

/**
 * This class overrides the solve method from Discretization with the Gauss Seidel method.
*/

class GaussSeidel : public PressureSolver
{
    public:
        //constructor
        GaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);
    
        //solve the Poisson problem for the pressure, using the Gauss-Seidel method
        void solve();
};

