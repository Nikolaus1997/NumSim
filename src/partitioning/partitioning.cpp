#include "partitioning.h"
#include <cmath>

Partitioning::Partitioning(std::array<int, 2> nCellsGlobal)
    : nCellsGlobal_(nCellsGlobal), Decomposition_(std::array<int, 2> ()) 
{
    //Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks_);

    //Get the ID of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo_);

    MPI_Dims_create(nRanks_, 2, Decomposition_.data());

    //Berechnung der lokalen Zellenanzahl in interger wird abgerundet
    nCellsLocal_[0] = nCellsGlobal_[0] / Decomposition_[0];
    nCellsLocal_[1] = nCellsGlobal_[1] / Decomposition_[1];

    //Berechnung der Position der Subdomain
    nodeposition_[0] = getColumnIndex(ownRankNo_);
    nodeposition_[1] = getRowIndex(ownRankNo_);   
    
    //Berechnung des Rests der Zellenanzahl
    int remainderCol = nCellsGlobal_[0] % Decomposition_[0];
    int remainderRow = nCellsGlobal_[1] % Decomposition_[1];

    //Berechnung der globalen Startzelle der Subdomain als faktor
    int columnDecomposition   = ownRankNo_% Decomposition_[0];
    int rowDecomposistion     = floor(ownRankNo_ / Decomposition_[0]);

    //Berechnung der globalen Startzelle
    int columnStart = columnDecomposition * nCellsLocal_[0];
    for(int i = 1; i<nodeposition_[0];i++)
    {
        if(i<= remainderCol) {
        columnStart++;
    }
    }

    //Berechnung der globalen Startzelle
    int rowStart = rowDecomposistion * nCellsLocal_[1];
    for(int i = 1; i<nodeposition_[1];i++)
    {
        if(i<= remainderRow) {
            rowStart++;
        }
    }
    //erweitern der lokalen Zellenanzahl um 1 falls nicht perfekte Aufteilung für die domains die am Anfang sind
    if(nodeposition_[0] <= remainderCol) {
        nCellsLocal_[0] += 1;
    }
    //erweitern der lokalen Zellenanzahl um 1 falls nicht perfekte Aufteilung für die domains die am Anfang sind
    if(nodeposition_[1] <= remainderRow) {
        nCellsLocal_[1] += 1;
    }

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



