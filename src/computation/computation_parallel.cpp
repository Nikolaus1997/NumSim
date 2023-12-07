#include "computation_parallel.h"
#include "partitioning/partitioning.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/sor.h"
#include "settings/settings.h"
#include "discretization/discretization.h"
#include "discretization/donor_cell.h"
#include "discretization/central_differences.h"
#include "pressure_solver/sor_red_black.h"

// initialize the computation object, parse the settings from file that is given as the only command line argument
void ComputationParallel::initialize(std::string filename)
{
    // load settings
    settings_ = Settings();
    settings_.loadFromFile(filename);
    meshWidth_[0] = settings_.physicalSize[0] / settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1] / settings_.nCells[1];

    // initialize partitioning
    partitioning_ = std::make_shared<Partitioning>(settings_.nCells);

    // initialize discretization
    if(partitioning_->ownRankNo()==0){
    std::cout << "===========================================================================================================================" << std::endl;
    }
    if (settings_.useDonorCell)
    {
        if(partitioning_->ownRankNo()==0)
            std::cout << "Using DonorCell..." << std::endl;
        
        discretization_ = std::make_shared<DonorCell>(partitioning_, meshWidth_, settings_.alpha);
    }
    else
    {
        if(partitioning_->ownRankNo()==0)
            std::cout << "Using CentralDifferences..." << std::endl;
        discretization_ = std::make_shared<CentralDifferences>(partitioning_, meshWidth_);
    }

    // initialize the pressure solver
    if (settings_.pressureSolver == "SOR")
    {
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon,
                                                settings_.maximumNumberOfIterations, settings_.omega);
    }
    else if (settings_.pressureSolver == "GaussSeidel")
    {
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon,
                                                        settings_.maximumNumberOfIterations);
    }else if(settings_.pressureSolver == "SORRedBlack")
    {   
        std::cout << "Using SORRedBlack..." << std::endl;
        pressureSolver_ = std::make_unique<SORRedBlack>(discretization_, settings_.epsilon,
                                                        settings_.maximumNumberOfIterations, settings_.omega, partitioning_ );
    }
    else
    {
        std::cout << "Solver not found!" << std::endl;
    }

    // initialize output writer
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discretization_,partitioning_);
    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, partitioning_);
}

// run the whole simulation until tend
void ComputationParallel::runSimulation()
{
    int t_i = 0;
    double time = 0.0;
    // main loop of the simulation over time
    while (time < settings_.endTime)
    {
        t_i++;
        applyBoundaryValues();
        applyBoundaryValuesFandG();
        dt_=0.05;
        computeTimeStepWidth();

        if (time + dt_ > settings_.endTime)
        {
            dt_ = settings_.endTime - time;
        }
        MPI_Allreduce(MPI_IN_PLACE, &dt_, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        time = time + dt_;
        computePreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        outputWriterText_->writeFile(time);
        outputWriterParaview_->writeFile(time);
    }
}

// compute the time step width dt from maximum velocities
void ComputationParallel::computeTimeStepWidth()
{
    double dt_diffusion = settings_.re / 2 / (1 / (discretization_->dx() * discretization_->dx()) + 1 / (discretization_->dy() * discretization_->dy()));

    double max_abs_u = 0.0;
    double max_abs_v = 0.0;
    double dt_u = 0.0;
    double dt_v = 0.0;

    // find maximum value of velocities u
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
        {
            if (fabs(discretization_->u(i, j)) > max_abs_u)
            {
                max_abs_u = fabs(discretization_->u(i, j));
            }
        }
    }

    // find maximum value of velocities v
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
    {
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
        {
            if (fabs(discretization_->v(i, j)) > max_abs_v)
            {
                max_abs_v = fabs(discretization_->v(i, j));
            }
        }
    }
    // set time step width
    MPI_Allreduce(MPI_IN_PLACE, &max_abs_u, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &max_abs_v, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dt_u = discretization_->dx() / max_abs_u;
    dt_v = discretization_->dy() / max_abs_v;
    dt_ = settings_.alpha * std::min(dt_u, dt_v);
    dt_ = std::min(dt_, dt_diffusion);
}

// set boundary values of u and v
void ComputationParallel::applyBoundaryValues()
{
    MPI_Request request;
    if(partitioning_->ownPartitionContainsTopBoundary())
    {
        // set boundary values u at top
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
                    discretization_->u(i, discretization_->uJEnd() - 1) =
                2.0 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 2);
        }
        // set boundary values v at top
        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
            discretization_->v(i, discretization_->vJEnd() - 1) = settings_.dirichletBcTop[1];
        }
    }else
    {
        partitioning_->mpiExchangeTop(discretization_->u(), request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        partitioning_->mpiExchangeTop(discretization_->v(), request); 
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    if(partitioning_->ownPartitionContainsBottomBoundary())
        {
                // set boundary values u at bottom
                for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
                    discretization_->u(i, discretization_->uJBegin()) =
                    2.0 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin()+1);
                }
                // set boundary values v at bottom
                for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
                    discretization_->v(i, discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
                }
        }else
        {
            partitioning_->mpiExchangeBottom(discretization_->u(), request);
            MPI_Wait(&request, MPI_STATUS_IGNORE);
            partitioning_->mpiExchangeBottom(discretization_->v(), request);
            MPI_Wait(&request, MPI_STATUS_IGNORE);
        }


    if(partitioning_->ownPartitionContainsLeftBoundary())
    {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
        {
            discretization_->u(discretization_->uIBegin(), j) = settings_.dirichletBcLeft[0];
        }

            // set boundary values v left 
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            discretization_->v(discretization_->vIBegin(), j) =
                    2.0 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin() + 1, j);
        }
    }else
    {
        partitioning_->mpiExchangeLeft(discretization_->u(), request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        partitioning_->mpiExchangeLeft(discretization_->v(), request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    if(partitioning_->ownPartitionContainsRightBoundary())
    {
        for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
        {
            discretization_->u(discretization_->uIEnd() - 1, j) = settings_.dirichletBcRight[0];
        }

        // set boundary values v left and right
        for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
            discretization_->v(discretization_->vIEnd() - 1, j) =
                    2.0 * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vIEnd() - 2, j);
        }
    }else
    {
        partitioning_->mpiExchangeRight(discretization_->u(), request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        partitioning_->mpiExchangeRight(discretization_->v(), request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    
    
}