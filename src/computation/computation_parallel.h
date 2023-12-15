#pragma once

#include "computation.h"
#include <array>
#include <memory>
#include <vector>
#include "output_writer/output_writer_paraview_parallel.h"
#include "output_writer/output_writer_text_parallel.h"



class ComputationParallel: public Computation{
    public:
        void initialize(std::string filename);
        void runSimulation();
    protected:
        void computeTimeStepWidth();
        void applyBoundaryValues();
        void setTestValues();
};