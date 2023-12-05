#pragma once

#include <memory>
#include <array>
#include <mpi.h>

class Partitioning {
    public:
        Partitioning(std::array<int, 2> nCellsGlobal_);
        int getRankID();
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
        bool checkLeftBoundary()                const;
        bool checkRightBoundary()               const;
        bool checkTopBoundary()                 const;
        bool checkBottomBoundary()              const;
        int calcRankID(int row, int column)     const;
        int getLeftAdjacentRankID()             const;
        int getRightAdjacentRankID()            const;
        int getTopAdjacentRankID()              const;
        int getBottomAdjacentRankID()           const;


    private:
        int nProcs_;
        int rankID_;
        std::array<int, 2> Decomposition_;
        std::array<int, 2> nCellsGlobal_;
        std::array<int, 2> nCellsLocal_;
        std::array<int, 2> nodestart_;
        std::array<int, 2> nodeposition_;
        int LeftNeighborRankID_;
        int RightNeighborRankID_;
        int TopNeighborRankID_;
        int BottomNeighborRankID_;

};