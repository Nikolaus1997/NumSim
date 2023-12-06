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
    //interior
    int pIBegin_in = discretization_->pIBegin();
    int pIEnd_in = discretization_->pIEnd()-1;
    int pJBegin_in = discretization_->pJBegin();
    int pJEnd_in = discretization_->pJEnd()-1;

    //calculate number of inner rows and columns
    int NColumns = (discretization_->pIEnd()-1)-(discretization_->pIBegin());
    int NRows = (discretization_->pJEnd()-1)-(discretization_->pJBegin());

    std::vector<double> TopRow(NRows, 0);
    std::vector<double> BottomRow(NRows, 0);
    std::vector<double> RightColumn(NColumns, 0);
    std::vector<double> LeftColumn(NColumns, 0);

    if (partitioning_->ownPartitionContainsTopBoundary()) {
        //set top boundary values
    } else {
        //send top row to top neighbour
        int End = discretization_->pIEnd()-1;
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd()-1; i++){
            TopRow.at(i - pJBegin_in) = discretization_->p(i,pJEnd_in);
        }
        MPI_Send(TopRow.data(), TopRow.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD);
        //receive ghost layer from top neighbour
        MPI_Recv(TopRow.data(), NRows, MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = pIBegin_in; i<pIEnd_in; i++){
            discretization_->p(i, discretization_->pJEnd()) = TopRow.at(i-pJBegin_in);
        }
    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        //set bottom boundary values
    } else {
        //send top row to top neighbour
        //receive ghost layer from bottom neighbour
    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        //set left boundary values
    } else {
        //send left column to left neighbour
        //receive ghost layer from left neighbour
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        //set right boundary values
    } else {
        //send right column to left neighbour
        //receive ghost layer from right neighbour
    }

}