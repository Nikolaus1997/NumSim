#include "pressure_solver/pressure_solver_parallel.h"
#include "partitioning/partitioning.h"
#include <mpi.h>

PressureSolverParallel::PressureSolverParallel(std::shared_ptr< Discretization >  discretization, 
                                                double epsilon, 
                                                int maximumNumberOfIterations, 
                                                std::shared_ptr<Partitioning> partitioning) : 
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), partitioning_(partitioning)
{
    N_= partitioning_->nCellsGlobal()[0]*partitioning_->nCellsGlobal()[1];
}
//TODO: implement communication
void PressureSolverParallel::communicateBoundaries(){
    //calculate number of inner rows and columns
    // int pIBegin_in = discretization_->pIBegin()+1;
    // int pIEnd_in = discretization_->pIEnd()-2;
    // int pJBegin_in = discretization_->pJBegin()+1;
    // int pJEnd_in = discretization_->pJEnd()-2;
    // int NColumns = pIEnd_in-pIBegin_in;
    // int NRows = pJEnd_in-pJBegin_in;

    MPI_Request request;

    // std::vector<double> TopRow(NRows, 0);
    // std::vector<double> BottomRow(NRows, 0);
    // std::vector<double> RightColumn(NColumns, 0);
    // std::vector<double> LeftColumn(NColumns, 0);

    if (partitioning_->ownPartitionContainsTopBoundary()) {
        //set top boundary values
        for(int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
        {
            //p top boundary condition
            discretization_->p(i,discretization_->pJEnd()-1)= discretization_->p(i,discretization_->pJEnd()-2);
        }
    } else {
        //send top row to top neighbour
        // for (int i = pIBegin_in; i <= pIEnd_in; i++){
        //     TopRow.at(i) = discretization_->p(i,pJEnd_in);
        // }
        // MPI_Send(TopRow.data(), TopRow.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD);
        // //receive ghost layer from top neighbour
        // MPI_Recv(TopRow.data(), NRows, MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // for (int i = pIBegin_in; i<=pIEnd_in; i++){
        //     discretization_->p(i, pJEnd_in+1) = TopRow.at(i);
        // }
        partitioning_->mpiExchangeTop(discretization_->p(), request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        //set bottom boundary values
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
        {
            //p bottom  boundary condition
            discretization_->p(i,discretization_->pJBegin()) = discretization_->p(i,discretization_->pJBegin()+1);
        }   
    } else {
        //send top row to top neighbour
        // for (int i = pIBegin_in; i <= pIEnd_in; i++){
        //     BottomRow.at(i) = discretization_->p(i,pJBegin_in);
        // }
        // MPI_Send(BottomRow.data(), BottomRow.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, MPI_COMM_WORLD);
        // //receive ghost layer from bottom neighbour
        // MPI_Recv(BottomRow.data(), NRows, MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // for (int i = pIBegin_in; i<=pIEnd_in; i++){
        //     discretization_->p(i, pIBegin_in-1) = BottomRow.at(i);
        // }
        partitioning_->mpiExchangeBottom(discretization_->p(), request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        //set left boundary values
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
        {
            //p left boundary condition
            discretization_->p(discretization_->pIBegin(),j) = discretization_->p(discretization_->pIBegin()+1, j);
        }
    } else {
        //send left column to left neighbour
        // for (int j = pJBegin_in; j <= pJEnd_in; j++){
        //     LeftColumn.at(j) = discretization_->p(pIBegin_in,j);
        // }
        // MPI_Send(LeftColumn.data(), LeftColumn.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, MPI_COMM_WORLD);
        // //receive ghost layer from left neighbour
        // MPI_Recv(LeftColumn.data(), NColumns, MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // for (int j = pJBegin_in; j<=pJEnd_in; j++){
        //     discretization_->p(pJBegin_in-1, j) = LeftColumn.at(j);
        // }
        partitioning_->mpiExchangeLeft(discretization_->p(), request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        //set right boundary values
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
        {
            //p right boundary condition
            discretization_->p(discretization_->pIEnd()-1, j)= discretization_->p(discretization_->pIEnd()-2, j);
        }
    } else {
        // //send right column to left neighbour
        // for (int j = pJBegin_in; j <= pJEnd_in; j++){
        //     RightColumn.at(j) = discretization_->p(pIEnd_in,j);
        // }
        // MPI_Send(RightColumn.data(), RightColumn.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, MPI_COMM_WORLD);
        // //receive ghost layer from right neighbour
        // MPI_Recv(LeftColumn.data(), NColumns, MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // for (int j = pJBegin_in; j<=pJEnd_in; j++){
        //     discretization_->p(pJBegin_in-1, j) = LeftColumn.at(j);
        // }
        partitioning_->mpiExchangeRight(discretization_->p(), request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

}

void PressureSolverParallel::computeResiduum()
{
    double dxdx = pow(discretization_->dx(),2);
    double dydy = pow(discretization_->dy(),2);
    residuum_ = 0.;
        for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
            for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
                    double d2pdxx = (discretization_->p(i-1,j) - 2 * discretization_->p(i, j) + discretization_->p(i+1,j))/dxdx;
                    double d2pdyy = (discretization_->p(i,j-1) - 2 * discretization_->p(i, j) + discretization_->p(i,j+1))/dydy; 

                    residuum_ += pow(d2pdxx + d2pdyy - discretization_->rhs(i,j),2);
            } 
            
        }
}