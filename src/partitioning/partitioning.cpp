#include "partitioning.h"
#include <cmath>
#include <iostream>

Partitioning::Partitioning(std::array<int, 2> nCellsGlobal)
    : nCellsGlobal_(nCellsGlobal), Decomposition_(std::array<int, 2> ()) 
{
    //Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks_);

    //Get the ID of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo_);

    MPI_Dims_create(nRanks_, 2, Decomposition_.data());
    std::cout<<"Decomposition: "<<Decomposition_[0]<<" "<<Decomposition_[1]<<std::endl;
    //Berechnung der lokalen Zellenanzahl in interger wird abgerundet
    nCellsLocal_[0] = nCellsGlobal_[0] / Decomposition_[0];
    nCellsLocal_[1] = nCellsGlobal_[1] / Decomposition_[1];

    //Berechnung der Position der Subdomain
    nodeposition_[0] = ownRankNo_ % Decomposition_[0]+1;
    nodeposition_[1] = (int) ownRankNo_ / Decomposition_[0]+1;   
    
    //Berechnung des Rests der Zellenanzahl
    int remainderCol = nCellsGlobal_[0] % Decomposition_[0];
    int remainderRow = nCellsGlobal_[1] % Decomposition_[1];

    int columnStart = 0;
    for (int i = 1; i < nodeposition_[0]; i++) {
        columnStart += nCellsLocal_[0];
        if (i <= remainderCol)
            columnStart++;
    }

    int rowStart = 0;
    for (int j = 1; j < nodeposition_[1]; j++) {
        rowStart += nCellsLocal_[1];
        if (j <= remainderRow)
            rowStart++;
    }
    if (nodeposition_[0] <= remainderCol)
        nCellsLocal_[0]++;
    if (nodeposition_[1] <= remainderRow)
        nCellsLocal_[1]++;

    std::cout<<"Rank: "<< ownRankNo_ <<" Größe: "<< nCellsLocal_[0] <<" | "<< nCellsLocal_[1]<<std::endl;

    //setzen der Nachbarn
    if( ownPartitionContainsLeftBoundary()) {
        LeftNeighborRankID_ = ownRankNo_;
    } else {
        LeftNeighborRankID_ = calcRankID(nodeposition_[0]-1, nodeposition_[1]);
    }
    if(ownPartitionContainsRightBoundary()) {
        RightNeighborRankID_ = ownRankNo_;
    } else {
        RightNeighborRankID_ = calcRankID(nodeposition_[0]+1, nodeposition_[1]);
    }
    if(ownPartitionContainsBottomBoundary()) {
        BottomNeighborRankID_ = ownRankNo_;
    } else {
        BottomNeighborRankID_ = calcRankID(nodeposition_[0], nodeposition_[1]-1);
    }
    if(ownPartitionContainsTopBoundary()) {
        TopNeighborRankID_ = ownRankNo_;
    } else {
        TopNeighborRankID_ = calcRankID(nodeposition_[0], nodeposition_[1]+1);
    }

    nodeOffset_ = {columnStart, rowStart};


    //MPI_Type_vector(nCellsLocal_[1], 1, nCellsLocal_[0], MPI_DOUBLE, &column);
    //MPI_Type_commit(&column);
    //MPI_Type_commit(&column);

    //MPI_Type_vector(nCellsLocal_[1], 1, 1, MPI_DOUBLE, &row);
    //MPI_Type_commit(&row);
    //MPI_Type_commit(&row);
};

int Partitioning::ownRankNo() const{
    return ownRankNo_;
}

int Partitioning::getNProcs() {
    return nRanks_;
}


std::array<int, 2> Partitioning::getDecomposition() const {
    return Decomposition_;
}   

std::array<int, 2> Partitioning::nCellsGlobal() const{
    return nCellsGlobal_;
}

std::array<int, 2> Partitioning::nCellsLocal() const{
    return nCellsLocal_;
}

int Partitioning::getColumnIndex(int rankID) const{
    return calcColumn(rankID)+1;
}

int Partitioning::calcColumn(int rankID) const {
    return rankID % Decomposition_[0];
}

int Partitioning::calcRow(int rankID) const {
    return floor(rankID / Decomposition_[0]);
}

int Partitioning::getRowIndex(int rankID) const{
    return calcRow(rankID)+1;
}

std::array<int, 2> Partitioning::nodeOffset() const{
    return nodeOffset_;
}

int Partitioning::getDecompositionColumnOrigin() const{
    return 1;
}

int Partitioning::getDecompositionRowOrigin()    const{
    return 1;
}

int Partitioning::getDecompositionColumnEnd()    const{
    return Decomposition_[0];
}

int Partitioning::getDecompositionRowEnd()       const{
    return Decomposition_[1];
}

bool Partitioning::ownPartitionContainsLeftBoundary()const{
    return  nodeposition_[0]== getDecompositionColumnOrigin();
}

bool Partitioning::ownPartitionContainsRightBoundary() const{
    return  nodeposition_[0]== getDecompositionColumnEnd();
}

int Partitioning::leftNeighbourRankNo() const
{
    return LeftNeighborRankID_;
}

int Partitioning::rightNeighbourRankNo() const
{
    return RightNeighborRankID_;
}

int Partitioning::topNeighbourRankNo() const
{
    return TopNeighborRankID_;
}

int Partitioning::bottomNeighbourRankNo() const
{
    return BottomNeighborRankID_;
}   


bool Partitioning::ownPartitionContainsTopBoundary() const{
    return  nodeposition_[1]== getDecompositionRowEnd();
}

bool Partitioning::ownPartitionContainsBottomBoundary() const{
    return  nodeposition_[1]== getDecompositionRowOrigin();
}

//nehme immer den Block in x richtung (Decomposition_[0]) und multipliziere ihn mit der anzahl von reihen
//addiere die Anzahl von Subdomain spalten die noch übgrig bleiben um bei der subdomain anzukommen 
//input sind hier die globalen subdomain indeces
int Partitioning::calcRankID(int column, int row) const{
    return (column-1) + (row-1) * Decomposition_[0];
}


/*void Partitioning::mpiExchangeAll(Array2D data, MPI_Request &request) const{
    mpiExchangeTop(data, request);
    mpiExchangeRight(data, request);
    mpiExchangeLeft(data, request);
    mpiExchangeBottom(data, request);
}


void Partitioning::mpiExchangeTop(Array2D data, MPI_Request &request) const
{   
    if (topNeighbourRankNo() != ownRankNo()){
        MPI_Isend(data.data()+nCellsLocal_[0], 1, row, topNeighbourRankNo(), 00, MPI_COMM_WORLD, &request);
    }
    if (bottomNeighbourRankNo() != ownRankNo()){
        MPI_Irecv(data.data()+(nCellsLocal_[1] - 1) * nCellsLocal_[0], 1, row, bottomNeighbourRankNo(), 00, MPI_COMM_WORLD, &request);
    }
}

void Partitioning::mpiExchangeRight(Array2D data, MPI_Request &request) const
{
    if(rightNeighbourRankNo() != ownRankNo()){
        MPI_Isend(data.data()+1, 1, column, rightNeighbourRankNo(), 11, MPI_COMM_WORLD, &request);
    }
    if(leftNeighbourRankNo() != ownRankNo()){
        MPI_Irecv(data.data()+nCellsLocal_[0] -1, 1, column, leftNeighbourRankNo(), 11, MPI_COMM_WORLD, &request);
    }
}

void Partitioning::mpiExchangeLeft(Array2D data, MPI_Request &request) const
{
    if(leftNeighbourRankNo() != ownRankNo()){
        MPI_Isend(data.data()+nCellsLocal_[0] -2, 1, column, leftNeighbourRankNo(), 22, MPI_COMM_WORLD, &request);
    }
    if(rightNeighbourRankNo() != ownRankNo()){
        MPI_Irecv(data.data(), 1, column, rightNeighbourRankNo(), 22, MPI_COMM_WORLD, &request);
    }
}

void Partitioning::mpiExchangeBottom(Array2D data, MPI_Request &request) const
{
    if (bottomNeighbourRankNo() != ownRankNo()){
        MPI_Isend(data.data()+(nCellsLocal_[1] - 2) * nCellsLocal_[0], 1, row, bottomNeighbourRankNo(), 33, MPI_COMM_WORLD, &request);
    }
    if (topNeighbourRankNo() != ownRankNo()){
        MPI_Irecv(data.data(), 1, row, topNeighbourRankNo(), 33, MPI_COMM_WORLD, &request);
    }
} */