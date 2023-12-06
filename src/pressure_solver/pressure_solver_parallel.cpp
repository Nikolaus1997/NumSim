#include "pressure_solver/pressure_solver_parallel.h"
#include "partitioning/partitioning.h"

PressureSolverParallel::PressureSolverParallel(std::shared_ptr< Discretization >  discretization, 
                                                double epsilon, 
                                                int maximumNumberOfIterations, 
                                                std::shared_ptr<Partitioning> partitioning) : 
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), partitioning_(partitioning)
{
}

void PressureSolverParallel::communicateBoundaries(){
    //calculate number of inner rows and columns
    int pIBegin_in = discretization_->pIBegin()+1;
    int pIEnd_in = discretization_->pIEnd()-2;
    int pJBegin_in = discretization_->pJBegin()+1;
    int pJEnd_in = discretization_->pJEnd()-2;
    int NColumns = pIEnd_in-pIBegin_in;
    int NRows = pJEnd_in-pJBegin_in;

    std::vector<double> TopRow(NRows, 0);
    std::vector<double> BottomRow(NRows, 0);
    std::vector<double> RightColumn(NColumns, 0);
    std::vector<double> LeftColumn(NColumns, 0);

    if (partitioning_->ownPartitionContainsTopBoundary()) {
        //set top boundary values
    } else {
        //send top row to top neighbour
        for (int i = pIBegin_in; i <= pIEnd_in; i++){
            TopRow.at(i) = discretization_->p(i,pJEnd_in);
        }
        MPI_Send(TopRow.data(), TopRow.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD);
        //receive ghost layer from top neighbour
        MPI_Recv(TopRow.data(), NRows, MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = pIBegin_in; i<=pIEnd_in; i++){
            discretization_->p(i, pJEnd_in+1) = TopRow.at(i);
        }
    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        //set bottom boundary values
    } else {
        //send top row to top neighbour
        for (int i = pIBegin_in; i <= pIEnd_in; i++){
            BottomRow.at(i) = discretization_->p(i,pJBegin_in);
        }
        MPI_Send(BottomRow.data(), BottomRow.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, MPI_COMM_WORLD);
        //receive ghost layer from bottom neighbour
        MPI_Recv(BottomRow.data(), NRows, MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = pIBegin_in; i<=pIEnd_in; i++){
            discretization_->p(i, pIBegin_in-1) = BottomRow.at(i);
        }
    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        //set left boundary values
    } else {
        //send left column to left neighbour
        for (int j = pJBegin_in; j <= pJEnd_in; j++){
            LeftColumn.at(j) = discretization_->p(pIBegin_in,j);
        }
        MPI_Send(LeftColumn.data(), LeftColumn.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, MPI_COMM_WORLD);
        //receive ghost layer from left neighbour
        MPI_Recv(LeftColumn.data(), NColumns, MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int j = pJBegin_in; j<=pJEnd_in; j++){
            discretization_->p(pJBegin_in-1, j) = LeftColumn.at(j);
        }
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        //set right boundary values
    } else {
        //send right column to left neighbour
        for (int j = pJBegin_in; j <= pJEnd_in; j++){
            RightColumn.at(j) = discretization_->p(pIEnd_in,j);
        }
        MPI_Send(RightColumn.data(), RightColumn.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, MPI_COMM_WORLD);
        //receive ghost layer from right neighbour
        MPI_Recv(LeftColumn.data(), NColumns, MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int j = pJBegin_in; j<=pJEnd_in; j++){
            discretization_->p(pJBegin_in-1, j) = LeftColumn.at(j);
        }
    }

}