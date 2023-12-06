#pragma once

#include "computation.h"
#include <array>
#include <memory>
#include "output_writer/output_writer_paraview_parallel.h"
#include "partitioning/partitioning.h"

class ComputationParallel: public Computation{
    public:
        void initialize(std::string filename);
        void runSimulation();
    protected:
        void computeTimeStepWidth();
        void applyBoundaryValues();

    std::shared_ptr<Partitioning> partitioning_;
};