#pragma once
#include "pressure_solver/pressure_solver_parallel.h"

class SORRedBlack : public PressureSolverParallel {
    public:
    //constructor
    SORRedBlack(std::shared_ptr<Discretization> discretization, 
                    double epsilon, 
                    int maximumNumberOfIterations, 
                    double omega,
                    std::shared_ptr<Partitioning> partitioning);
    
    //solve the poisson problem using the red black algorithm
    void solve();

    private:
        double omega_;
};