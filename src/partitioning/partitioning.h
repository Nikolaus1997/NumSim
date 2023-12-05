#pragma once
#include <array>
#include <memory>
#include <array>
#include <mpi.h>

class Partitioning {
    public:
        Partitioning(std::array<int, 2> nCellsGlobal_);
        int ownRankNo() const;
        int getNProcs();
        std::array<int, 2> getDecomposition()   const;
        std::array<int, 2> getNCellsGlobal()    const;
        std::array<int, 2> getNCellsLocal()     const;
        int calcColumn      (int rankID)        const;
        int getColumnIndex  (int rankID)        const;
        int calcRow         (int rankID)        const;
        int getRowIndex     (int rankID)        const;
        std::array<int, 2> getNodeIndeces()     const;
        int getDecompositionColumnOrigin()      const;                     
        int getDecompositionRowOrigin()         const;   
        int getDecompositionColumnEnd()         const;
        int getDecompositionRowEnd()            const;   
        //! if the own partition has part of the bottom boundary of the whole domain
        bool ownPartitionContainsBottomBoundary() const;
        //! if the own partition has part of the top boundary of the whole domain
        //! used in OutputWriterParaviewParallel
        bool ownPartitionContainsTopBoundary() const;
        //! if the own partition has part of the left boundary of the whole domain
        bool ownPartitionContainsLeftBoundary() const;
        //! if the own partition has part of the right boundary of the whole domain
        //! used in OutputWriterParaviewParallel
        bool ownPartitionContainsRightBoundary() const;
        //! get the rank no of the left neighbouring rank
        int leftNeighbourRankNo() const;
        //! get the rank no of the right neighbouring rank
        int rightNeighbourRankNo() const;
        //! get the rank no of the top neighbouring rank
        int topNeighbourRankNo() const;
        //! get the rank no of the bottom neighbouring rank
        int bottomNeighbourRankNo() const;
        int calcRankID(int row, int column)     const;
        std::array<int,2> nodeOffset() const;


    private:
        int nRanks_;
        int ownRankNo_;
        std::array<int, 2> Decomposition_;
        std::array<int, 2> nCellsGlobal_;
        std::array<int, 2> nCellsLocal_;
        std::array<int, 2> nodeOffset_;
        std::array<int, 2> nodeposition_;
        int LeftNeighborRankID_ ;
        int RightNeighborRankID_;
        int TopNeighborRankID_;
        int BottomNeighborRankID_;

};