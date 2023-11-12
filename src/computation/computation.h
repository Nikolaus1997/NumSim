#pragma once

#include "settings/settings.h"
#include "discretization/discretization.h"
#include "discretization/donor_cell.h"
#include "discretization/central_differences.h"
#include "pressure_solver/pressure_solver.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"

#include <memory>
#include <cmath>
#include <algorithm>


class Computation
{
    public:
        void initialize(int argc, char *argv[]);

        void runSimulation();
    
    private:
        void computeTimeStepWidth();

        void applyBoundaryValues();

        void computePreliminaryVelocities();

        void computeRightHandSide();

        void computePressure();

        void computeVelocities();


        Settings settings_;
        std::shared_ptr<Discretization> discretization_;
        std::unique_ptr<PressureSolver> pressureSolver_;
        std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
        std::unique_ptr<OutputWriterText> outputWriterText_;
        std::array<double, 2> meshWidth_;
        double dt_;    
};