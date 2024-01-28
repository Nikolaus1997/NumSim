#include "discretization/lbm_discretization.h"

LbmDiscretization::LbmDiscretization(std::array<int, 2> nCells, std::array<double, 2> meshWidth):
                                    CentralGrid(nCells, meshWidth)
{
};

void LbmDiscretization::calcRho(int i, int j)
{
    double sum = 0.0;
        for(int k; k < 9 ; k++)
            {
                sum = sum + pdf(i,j,k);
            }
    rho(i,j) = sum;
};