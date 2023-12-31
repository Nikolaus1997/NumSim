#pragma once
#include "pressure_solver/pressure_solver.h"
#include "discretization/discretization.h"
#include "partitioning/partitioning.h"

class PressureSolverParallel : public PressureSolver {
public:
    //constructor
    PressureSolverParallel (std::shared_ptr<Discretization> discretization, 
                    double epsilon, 
                    int maximumNumberOfIterations, 
                    std::shared_ptr<Partitioning> partitioning);
    virtual void computeResiduum() ;

    protected:
    //communicate pressure boundaries
    void communicateBoundaries();

    std::shared_ptr<Partitioning> partitioning_;

};