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

/**
 * This class contains the main loop over all time steps of the simulation and all methods that are called in this loop.
*/

class Computation
{
    public:
        //initialize the computation object, parse the settings from file that is given as the only command line argument 
        void initialize(std::string filename);

        //run the whole simulation until tend 
        void runSimulation();
    
    private:
        //compute the time step width dt from maximum velocities 
        void computeTimeStepWidth();

        //set boundary values of u and v to correct values 
        void applyBoundaryValues();

        //set boundary values of F and G to correct values 
        void applyBoundaryValuesFandG();

        //compute the preliminary velocities, F and G 
        void computePreliminaryVelocities();

        //compute the right hand side for the poisson problem for pressure
        void computeRightHandSide();

        //solve the Poisson equation for the pressure 
        void computePressure();

        //compute the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p 
        void computeVelocities();


        Settings settings_;
        std::shared_ptr<Discretization> discretization_;
        std::unique_ptr<PressureSolver> pressureSolver_;
        std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
        std::unique_ptr<OutputWriterText> outputWriterText_;
        std::array<double, 2> meshWidth_;
        double dt_;    
};