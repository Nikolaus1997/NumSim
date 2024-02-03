#pragma once

#include "discretization/central_grid.h"
#include <array>


class LbmDiscretization 
                : public CentralGrid
{
    public:
        LbmDiscretization(std::array<int, 3> nCells,
                  std::array<double, 3> meshWidth); 
        
        void calcRho(int i, int j, int k);
};