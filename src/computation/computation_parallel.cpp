#include "computation_parallel.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/sor.h"
#include "settings/settings.h"
#include "discretization/discretization.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"
#include "discretization/donor_cell.h"
#include "discretization/central_differences.h"

//initialize the computation object, parse the settings from file that is given as the only command line argument 
void Computation_Parallel::initialize(std::string filename)
{   
    //load settings
    settings_ = Settings();
    settings_.loadFromFile(filename);
    meshWidth_[0] = settings_.physicalSize[0]/settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1]/settings_.nCells[1];
    
    //initialize discretization
    std::cout<<"==========================================================================================================================="<<std::endl;
    if (settings_.useDonorCell) {
        std::cout<<"Using DonorCell..."<<std::endl;
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    }
    else {
        std::cout<<"Using CentralDifferences..."<<std::endl;
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    }
    
    //initialize the pressure solver
    if (settings_.pressureSolver == "SOR") {
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon,
                                                settings_.maximumNumberOfIterations, settings_.omega);
    }
    else if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon,
                                                settings_.maximumNumberOfIterations);

    } else {
        std::cout << "Solver not found!" << std::endl;
    }
    //initialize output writer
    //outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
    //outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
}

//run the whole simulation until tend 
void Computation_Parallel::runSimulation()
{
    int t_i = 0;
    double time = 0.0;
    //main loop of the simulation over time
    while(time<settings_.endTime){
        t_i++;
        applyBoundaryValues();
        applyBoundaryValuesFandG();;
        computeTimeStepWidth();

        if (time + dt_>settings_.endTime){
            dt_ = settings_.endTime-time;
        }
        time = time + dt_;
        computePreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        //outputWriterText_->writeFile(time);
        outputWriterParaview_->writeFile(time);
    }
}

//compute the time step width dt from maximum velocities 
void Computation::computeTimeStepWidth()
{
    double dt_diffusion = settings_.re / 2 / (1 / (discretization_->dx() * discretization_->dx()) + 1 / (discretization_->dy() * discretization_->dy()) );
    
    double max_abs_u = 0.0;
    double max_abs_v = 0.0;
    double dt_u = 0.0;
    double dt_v = 0.0;

    //find maximum value of velocities u
    for(int i = discretization_->uIBegin(); i <discretization_->uIEnd() ;i++){
        for(int j = discretization_->uJBegin(); j < discretization_->uJEnd() ;j++){
            if(fabs(discretization_->u(i,j)) > max_abs_u){
                max_abs_u = fabs(discretization_->u(i,j));
            }
        }
    }

    //find maximum value of velocities v
    for(int i = discretization_->vIBegin(); i <discretization_->vIEnd() ;i++){
        for(int j = discretization_->vJBegin(); j < discretization_->vJEnd() ;j++){
            if(fabs(discretization_->v(i,j)) > max_abs_v){
                max_abs_v = fabs(discretization_->v(i,j));
            }
        }
    }   
    //set time step width
    dt_u = discretization_->dx()/max_abs_u;
    dt_v = discretization_->dy()/max_abs_v;
    dt_ = 0.05;
 }

//set boundary values of u and v 
void Computation_Parallel::applyBoundaryValues()
{
   
   
}




