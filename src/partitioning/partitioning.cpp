#include "partitioning.h"
#include <cmath>

Partitioning::Partitioning(std::array<int, 2> nCellsGlobal)
    : nCellsGlobal_(nCellsGlobal), Decomposition_(std::array<int, 2> ()) 
{
    //Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs_);

    //Get the ID of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &rankID_);

    MPI_Dims_create(nProcs_, 2, Decomposition_.data());

    nCellsLocal_[0] = nCellsGlobal_[0] / Decomposition_[0];
    nCellsLocal_[1] = nCellsGlobal_[1] / Decomposition_[1];

    int remainderCol = nCellsGlobal_[0] % Decomposition_[0];
    int remainderRow = nCellsGlobal_[1] % Decomposition_[1];

    int columnDecomposition   = rankID_ % Decomposition_[0] ;
    int rowDecomposistion     = floor(rankID_ / Decomposition_[0]) ;

    //Berechnung der globalen Startzelle
    int columnStart = columnDecomposition * nCellsLocal_[0];
    columnStart += std::min(columnDecomposition, remainderCol);

    //Berechnung der globalen Startzelle
    int rowStart = rowDecomposistion * nCellsLocal_[1];
    rowStart += std::min(rowDecomposistion, remainderRow);
    
    //erweitern der lokalen Zellenanzahl um 1 falls nicht perfekte Aufteilung
    if(columnDecomposition < remainderCol) {
        nCellsLocal_[0] += 1;
    }
    //erweitern der lokalen Zellenanzahl um 1 falls nicht perfekte Aufteilung
    if(rowDecomposistion < remainderRow) {
        nCellsLocal_[1] += 1;
    }


};
