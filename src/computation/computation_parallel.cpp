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
    MPI_Request uBottomSendRequest;
    MPI_Request uTopSendRequest;
    MPI_Request uLeftSendRequest;
    MPI_Request uRightSendRequest;

    MPI_Request uBottomReceiveRequest;
    MPI_Request uTopReceiveRequest;
    MPI_Request uLeftReceiveRequest;
    MPI_Request uRightReceiveRequest;

    MPI_Request vBottomSendRequest;
    MPI_Request vTopSendRequest;
    MPI_Request vLeftSendRequest;
    MPI_Request vRightSendRequest;

    MPI_Request vBottomReceiveRequest;
    MPI_Request vTopReceiveRequest;
    MPI_Request vLeftReceiveRequest;
    MPI_Request vRightReceiveRequest;

    int uIBeginInter = discretization_->uIBegin()+1;
    int uIEndInter = discretization_->uIEnd()-2;    
    int uColSizeInter = uIEndInter-uIBeginInter;

    int uJBeginInter = discretization_->uJBegin()+1;
    int uJEndInter = discretization_->uJEnd()-2;
    int uRowSizeInter = uJEndInter-uJBeginInter;

    int vIBeginInter = discretization_->vIBegin()+1;
    int vIEndInter = discretization_->vIEnd()-2;    
    int vColSizeInter = vIEndInter-vIBeginInter;

    int vJBeginInter = discretization_->vJBegin()+1;
    int vJEndInter = discretization_->vJEnd()-2;
    int vRowSizeInter = uJEndInter-uJBeginInter;


    std::vector<double> uBottomSendBuffer(uColSizeInter,0);
    std::vector<double> uTopSendBuffer(uColSizeInter,0);
    std::vector<double> uLeftSendBuffer(uRowSizeInter,0);
    std::vector<double> uRightSendBuffer(uRowSizeInter,0);

    std::vector<double> uBottomReceiveBuffer(uColSizeInter,0);
    std::vector<double> uTopReceiveBuffer(uColSizeInter,0);
    std::vector<double> uLeftReceiveBuffer(uRowSizeInter,0);
    std::vector<double> uRightReceiveBuffer(uRowSizeInter,0);

    std::vector<double> vBottomSendBuffer(uColSizeInter,0);
    std::vector<double> vTopSendBuffer(uColSizeInter,0);
    std::vector<double> vLeftSendBuffer(uRowSizeInter,0);
    std::vector<double> vRightSendBuffer(uRowSizeInter,0);

    std::vector<double> vBottomReceiveBuffer(uColSizeInter,0);
    std::vector<double> vTopReceiveBuffer(uColSizeInter,0);
    std::vector<double> vLeftReceiveBuffer(uRowSizeInter,0);
    std::vector<double> vRightReceiveBuffer(uRowSizeInter,0);


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
        //write top row to buffer
        for (int i=uIBeginInter; i<uIEndInter; i++){
            uTopSendBuffer.at(i - uIBeginInter) = discretization_->u(i,uJEndInter);
        }
        //send top buffer to top neighbour
        MPI_Isend(uTopSendBuffer.data(), uTopSendBuffer.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 00, MPI_COMM_WORLD, &uTopSendRequest);
        //receive bottom row from top neighbour
        MPI_Irecv(uTopReceiveBuffer.data(), uTopReceiveBuffer.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 00, MPI_COMM_WORLD, &uTopReceiveRequest);
        //wait for send and receive calls to complete
        MPI_Wait(&uTopSendRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&uTopReceiveRequest, MPI_STATUS_IGNORE);
        //write bottom row from top neighbour into ghost layer on top
        for (int i=uIBeginInter; i<uIEndInter; i++){
            discretization_->u(i,discretization_->uJEnd()-1) = uTopReceiveBuffer.at(i - uIBeginInter); 
        }

        //write top row to buffer
        for (int i=vIBeginInter; i<vIEndInter; i++){
            vTopSendBuffer.at(i - vIBeginInter) = discretization_->v(i,vJEndInter);
        }
        //send top buffer to top neighbour
        MPI_Isend(vTopSendBuffer.data(), vTopSendBuffer.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 00, MPI_COMM_WORLD, &vTopSendRequest);
        //receive bottom row from top neighbour
        MPI_Irecv(vTopReceiveBuffer.data(), vTopReceiveBuffer.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 00, MPI_COMM_WORLD, &vTopReceiveRequest);
        //wait for send and receive calls to complete
        MPI_Wait(&vTopSendRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&vTopReceiveRequest, MPI_STATUS_IGNORE);
        //write bottom row from top neighbour into ghost layer on top
        for (int i=vIBeginInter; i<vIEndInter; i++){
            discretization_->v(i,discretization_->vJEnd()-1) = vTopReceiveBuffer.at(i - vIBeginInter); 
        }

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
        //write bottom row to buffer
        for (int i=uIBeginInter; i<uIEndInter; i++){
            uTopSendBuffer.at(i - uIBeginInter) = discretization_->u(i,uJBeginInter);
        }
        //send bottom buffer to bottom neighbour
        MPI_Isend(uBottomSendBuffer.data(), uBottomSendBuffer.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 00, MPI_COMM_WORLD, &uBottomSendRequest);
        //receive bottom row from top neighbour
        MPI_Irecv(uBottomReceiveBuffer.data(), uBottomReceiveBuffer.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 00, MPI_COMM_WORLD, &uBottomReceiveRequest);
        //wait for send and receive calls to complete
        MPI_Wait(&uBottomSendRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&uBottomReceiveRequest, MPI_STATUS_IGNORE);
        //write bottom row from top neighbour into ghost layer on top
        for (int i=uIBeginInter; i<uIEndInter; i++){
            discretization_->u(i,discretization_->uJBegin()) = uBottomReceiveBuffer.at(i - uIBeginInter); 
        }

        //write top row to buffer
        for (int i=vIBeginInter; i<vIEndInter; i++){
            vBottomSendBuffer.at(i - vIBeginInter) = discretization_->v(i,vJEndInter);
        }
        //send top buffer to top neighbour
        MPI_Isend(vBottomSendBuffer.data(), vBottomSendBuffer.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 00, MPI_COMM_WORLD, &vBottomSendRequest);
        //receive bottom row from top neighbour
        MPI_Irecv(vBottomReceiveBuffer.data(), vBottomReceiveBuffer.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 00, MPI_COMM_WORLD, &vBottomReceiveRequest);
        //wait for send and receive calls to complete
        MPI_Wait(&vBottomSendRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&vBottomReceiveRequest, MPI_STATUS_IGNORE);
        //write bottom row from top neighbour into ghost layer on top
        for (int i=vIBeginInter; i<vIEndInter; i++){
            discretization_->v(i,discretization_->vJBegin()) = vTopReceiveBuffer.at(i - vIBeginInter); 
        }

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
        //write left row to buffer
        for (int j=uJBeginInter; j<uJEndInter; j++){
            uLeftSendBuffer.at(j - uJBeginInter) = discretization_->u(uIBeginInter,j);
        }
        //send left buffer to left neighbour
        MPI_Isend(uLeftSendBuffer.data(), uLeftSendBuffer.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 00, MPI_COMM_WORLD, &uLeftSendRequest);
        //receive right column from left neighbour
        MPI_Irecv(uLeftReceiveBuffer.data(), uLeftReceiveBuffer.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 00, MPI_COMM_WORLD, &uLeftReceiveRequest);
        //wait for send and receive calls to complete
        MPI_Wait(&uLeftSendRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&uLeftReceiveRequest, MPI_STATUS_IGNORE);
        //write right column from left neighbour into ghost layer on left
        for (int j=uJBeginInter; j<uJEndInter; j++){
            discretization_->u(discretization_->uIBegin(),j) = uLeftReceiveBuffer.at(j - uJBeginInter); 
        }

        //write left row to buffer
        for (int j=vJBeginInter; j<vJEndInter; j++){
            vLeftSendBuffer.at(j - vJBeginInter) = discretization_->v(vIBeginInter,j);
        }
        //send left buffer to left neighbour
        MPI_Isend(vLeftSendBuffer.data(), vLeftSendBuffer.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 00, MPI_COMM_WORLD, &vLeftSendRequest);
        //receive right column from left neighbour
        MPI_Irecv(vLeftReceiveBuffer.data(), vLeftReceiveBuffer.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 00, MPI_COMM_WORLD, &vLeftReceiveRequest);
        //wait for send and receive calls to complete
        MPI_Wait(&vLeftSendRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&vLeftReceiveRequest, MPI_STATUS_IGNORE);
        //write right column from left neighbour into ghost layer on left
        for (int j=vJBeginInter; j<vJEndInter; j++){
            discretization_->v(discretization_->vIBegin(),j) = vLeftReceiveBuffer.at(j - vJBeginInter); 
        }

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
        //write right column to buffer
        for (int j=uJBeginInter; j<uJEndInter; j++){
            uRightSendBuffer.at(j - uJBeginInter) = discretization_->u(uIEndInter,j);
        }
        //send right buffer to right neighbour
        MPI_Isend(uRightSendBuffer.data(), uRightSendBuffer.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 00, MPI_COMM_WORLD, &uRightSendRequest);
        //receive left column from right neighbour
        MPI_Irecv(uRightReceiveBuffer.data(), uRightReceiveBuffer.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 00, MPI_COMM_WORLD, &uRightReceiveRequest);
        //wait for send and receive calls to complete
        MPI_Wait(&uRightSendRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&uRightReceiveRequest, MPI_STATUS_IGNORE);
        //write left column from right neighbour into ghost layer on right
        for (int j=uJBeginInter; j<uJEndInter; j++){
            discretization_->u(discretization_->uIEnd()-1,j) = uRightReceiveBuffer.at(j - uJBeginInter); 
        }
        //write right column to buffer
        for (int j=vJBeginInter; j<vJEndInter; j++){
            vRightSendBuffer.at(j - vJBeginInter) = discretization_->v(uIEndInter,j);
        }
        //send right buffer to right neighbour
        MPI_Isend(vRightSendBuffer.data(), vRightSendBuffer.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 00, MPI_COMM_WORLD, &vRightSendRequest);
        //receive left column from right neighbour
        MPI_Irecv(vRightReceiveBuffer.data(), vRightReceiveBuffer.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 00, MPI_COMM_WORLD, &vRightReceiveRequest);
        //wait for send and receive calls to complete
        MPI_Wait(&vRightSendRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&vRightReceiveRequest, MPI_STATUS_IGNORE);
        //write left column from right neighbour into ghost layer on right
        for (int j=vJBeginInter; j<vJEndInter; j++){
            discretization_->v(discretization_->vIEnd()-1,j) = uRightReceiveBuffer.at(j - vJBeginInter); 
        }
    }


    }

    
    
