#pragma once
#include "pressure_solver.h"
#include <cmath>
#include <algorithm>

class GaussSeidel : public PressureSolver
{
    public:
        GaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);
    
        void solve();
};

