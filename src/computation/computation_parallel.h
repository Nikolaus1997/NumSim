#pragma once

#include "computation.h"


#include <array>
#include <memory>


class Computation_Parallel: public Computation{
    public:
        virtual void initialize(std::string filename);
        virtual void runSimulation();
    protected:
        void computeTimeStepWidth();
        void applyBoundaryValues();
};