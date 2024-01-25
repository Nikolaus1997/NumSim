#pragma once

#include "discretization/central_grid.h"
#include <array>

class LbmDiscretization 
                : CentralGrid
{
    public:
        LbmDiscretization(std::array<int, 2> nCells,
                  std::array<double, 2> meshWidth);

        
};