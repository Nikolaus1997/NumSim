#include "pressure_solver/pressure_solver_parallel.h"
#include "partitioning/partitioning.h"
#include <mpi.h>
#include <iostream>

PressureSolverParallel::PressureSolverParallel(std::shared_ptr< Discretization >  discretization, 
                                                double epsilon, 
                                                int maximumNumberOfIterations, 
                                                std::shared_ptr<Partitioning> partitioning) : 
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), partitioning_(partitioning)
{

}
//TODO: implement communication
void PressureSolverParallel::communicateBoundaries(){
    //calculate number of inner rows and columns
    int pIBegin_in = discretization_->pIBegin() + 1;
    int pIEnd_in = discretization_->pIEnd() - 1;
    int pJBegin_in = discretization_->pJBegin() + 1;
    int pJEnd_in = discretization_->pJEnd() - 1;
    int NColumns = pIEnd_in-pIBegin_in;
    int NRows = pJEnd_in-pJBegin_in;
    //std::cout << "NColumns " << NColumns << "NRows " << NRows << "pIBegin_in " << pIBegin_in << "pIEnd_in " << pIEnd_in << std::endl;

    MPI_Request topSendRequest;
    MPI_Request topReceiveRequest;
    MPI_Request bottomSendRequest;
    MPI_Request bottomReceiveRequest;
    MPI_Request leftSendRequest;
    MPI_Request leftReceiveRequest;
    MPI_Request rightSendRequest;
    MPI_Request rightReceiveRequest;

    //initialize buffer vectors
    std::vector<double> topSendBuffer(NColumns, 7);
    std::vector<double> topReceiveBuffer(NColumns, 7);
    std::vector<double> bottomSendBuffer(NColumns, 7);
    std::vector<double> bottomReceiveBuffer(NColumns, 7);
    std::vector<double> rightSendBuffer(NRows, 7);
    std::vector<double> rightReceiveBuffer(NRows, 7);
    std::vector<double> leftSendBuffer(NRows, 7);
    std::vector<double> leftReceiveBuffer(NRows, 7);

    if (partitioning_->ownPartitionContainsTopBoundary()) {
        //set top boundary values
        for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
        {
            //p top boundary condition
            discretization_->p(i,discretization_->pJEnd()-1)= discretization_->p(i,discretization_->pJEnd()-2);
        }
    } else {
        //write top row to buffer
        for (int i=pIBegin_in; i<pIEnd_in; i++){
            //std::cout <<partitioning_->ownRankNo()<<"  "<<topReceiveBuffer.at(i-pIBegin_in)<<std::endl;
            topSendBuffer.at(i) = discretization_->p(i,discretization_->pJEnd()-2);
        }
        //send top buffer to top neighbour
        MPI_Isend(topSendBuffer.data(), topSendBuffer.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 1, MPI_COMM_WORLD, &topSendRequest);
        //receive bottom row from top neighbour
        MPI_Irecv(topReceiveBuffer.data(), topReceiveBuffer.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 1, MPI_COMM_WORLD, &topReceiveRequest);

    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        //set bottom boundary values
        for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
        {
            //p bottom  boundary condition
            discretization_->p(i,discretization_->pJBegin()) = discretization_->p(i,discretization_->pJBegin()+1);
        }   
    } else {
        //write bottom row to buffer
        for (int i=pIBegin_in; i<pIEnd_in; i++){
            bottomSendBuffer.at(i) = discretization_->p(i,discretization_->pJBegin()+1);
        }
        //send bottom buffer to bottom neighbour
        MPI_Isend(bottomSendBuffer.data(), bottomSendBuffer.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 1, MPI_COMM_WORLD, &bottomSendRequest);
        //receive top row from bottom neighbour
        MPI_Irecv(bottomReceiveBuffer.data(), bottomReceiveBuffer.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 1, MPI_COMM_WORLD, &bottomReceiveRequest);

    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        //set left boundary values
        for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
        {
            //p left boundary condition
            discretization_->p(discretization_->pIBegin(),j) = discretization_->p(discretization_->pIBegin()+1, j);
        }
    } else {
        //write left row to buffer
        for (int j=pJBegin_in; j<pJEnd_in; j++){
            leftSendBuffer.at(j) = discretization_->p(discretization_->pIBegin()+1,j);
        }
        //send left buffer to left neighbour
        MPI_Isend(leftSendBuffer.data(), leftSendBuffer.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 1, MPI_COMM_WORLD, &leftSendRequest);
        //receive right column from left neighbour
        MPI_Irecv(leftReceiveBuffer.data(), leftReceiveBuffer.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(),1, MPI_COMM_WORLD, &leftReceiveRequest);

    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        //set right boundary values
        for(int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
        {
            //p right boundary condition
            discretization_->p(discretization_->pIEnd()-1, j)= discretization_->p(discretization_->pIEnd()-2, j);
        }
    } else {
        //write right column to buffer
        for (int j=pJBegin_in; j<pJEnd_in; j++){
            rightSendBuffer.at(j) = discretization_->p(discretization_->pIEnd()-2,j);
        }
        //send right buffer to right neighbour
        MPI_Isend(rightSendBuffer.data(), rightSendBuffer.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 1, MPI_COMM_WORLD, &rightSendRequest);
        //receive left column from right neighbour
        MPI_Irecv(rightReceiveBuffer.data(), rightReceiveBuffer.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 1, MPI_COMM_WORLD, &rightReceiveRequest);

    }

    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        //wait for send and receive calls to complete
        MPI_Wait(&topSendRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&topReceiveRequest, MPI_STATUS_IGNORE);
        //write bottom row from top neighbour into ghost layer on top
        for (int i=pIBegin_in; i<pIEnd_in; i++){
            //std::cout << topReceiveBuffer.at(i) << std::endl;
            discretization_->p(i,discretization_->pJEnd()-1) = topReceiveBuffer.at(i); 
        }
    } 

    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        //wait for send and receive calls to complete
        MPI_Wait(&bottomSendRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&bottomReceiveRequest, MPI_STATUS_IGNORE);
        //write bottom row from top neighbour into ghost layer on top
        for (int i=pIBegin_in; i<pIEnd_in; i++){
            //std::cout << bottomReceiveBuffer.at(i) << std::endl;
            discretization_->p(i,discretization_->pJBegin()) = bottomReceiveBuffer.at(i); 
        }
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        //wait for send and receive calls to complete
        MPI_Wait(&leftSendRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&leftReceiveRequest, MPI_STATUS_IGNORE);
        //write right column from left neighbour into ghost layer on left
        for (int j=pJBegin_in; j<pJEnd_in; j++){
            //std::cout << leftReceiveBuffer.at(j) << std::endl;
            discretization_->p(discretization_->pIBegin(),j) = leftReceiveBuffer.at(j); 
        }
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        //wait for send and receive calls to complete
        MPI_Wait(&rightSendRequest, MPI_STATUS_IGNORE);
        MPI_Wait(&rightReceiveRequest, MPI_STATUS_IGNORE);
        //write left column from right neighbour into ghost layer on right
        for (int j=pJBegin_in; j<pJEnd_in; j++){
            //std::cout << rightReceiveBuffer.at(j) << std::endl;
            discretization_->p(discretization_->pIEnd()-1,j) = rightReceiveBuffer.at(j); 
        }
    }
}

void PressureSolverParallel::computeResiduum()
{
    N_= partitioning_->nCellsGlobal()[0]*partitioning_->nCellsGlobal()[1];
    double holder = 0.0;
    double dxdx = pow(discretization_->dx(),2);
    double dydy = pow(discretization_->dy(),2);
        for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
            for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
                    double d2pdxx = (discretization_->p(i-1,j) - 2 * discretization_->p(i, j) + discretization_->p(i+1,j))/dxdx;
                    double d2pdyy = (discretization_->p(i,j-1) - 2 * discretization_->p(i, j) + discretization_->p(i,j+1))/dydy; 

                    residuum_ += pow(d2pdxx + d2pdyy - discretization_->rhs(i,j),2);
            } 
        }
    MPI_Allreduce(&residuum_, &holder,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    residuum_ = holder/N_;
}   