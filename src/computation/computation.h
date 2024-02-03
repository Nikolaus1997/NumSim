#pragma once
#include <memory>
#include <algorithm>
#include <cmath>
#include "settings/settings.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"
#include "discretization/lbm_discretization.h"
#include "storage/array2d.h"

/** This class handles the main simulation.
* It implements the time stepping scheme, computes all the terms and calls the pressure solver.
*/
class Computation {
public:
    /**
     * Initialize the computation object
     * 
     * Parse the settings from the parameter file that is given as the command line argument
     * It implements the time stepping scheme, computes all the terms and calls the pressure solver.
     */
    void initialize(std::string filename);

    /**
     * Run the whole simulation until tend
     */
    void runSimulation();

private:
    /**
     * Set the boundary values of the velocities (u, v)
     * 
     * Left and right boundaries should overwrite bottom and top boundaries
     */
   

    void LBMapplyBoundaryValues();

    void FillVelocitiesAndPressure();

    void Collision();
    
    void LBMBounceBack();

    void Streaming();

    Settings settings_;
    std::shared_ptr<LbmDiscretization> cdiscretization_;
    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
    std::unique_ptr<OutputWriterText> outputWriterText_;
    std::array<double, 3> meshWidth_;
    double dt_;
    double cs_;
    bool doLBM;
    double nu_;
    double tau_;
    double    tau_m ;
    double    tau_om_m ;
    double    tau_om_p ;
    double    tau_om   ;
};