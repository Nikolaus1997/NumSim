#pragma once
#include "pressure_solver/pressure_solver.h"

class SOR : public PressureSolver
{
    public:
        SOR (std::shared_ptr< Discretization > discretization, double epsilon, int maximumNumberOfIterations, double omega);
    
        void solve();
    
    private:
        double omega_;
}

