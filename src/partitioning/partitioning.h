#pragma once

#include <memory>
#include <array>
#include <mpi.h>

class Partitioning {
    public:
        Partitioning(std::array<int, 2> nCellsGlobal_);
    private:
        int nProcs_;
        int rankID_;
        std::array<int, 2> Decomposition_;
        std::array<int, 2> nCellsGlobal_;
        std::array<int, 2> nCellsLocal_;



};