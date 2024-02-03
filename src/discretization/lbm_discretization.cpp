#include "discretization/lbm_discretization.h"

LbmDiscretization::LbmDiscretization(std::array<int, 3> nCells, std::array<double, 3> meshWidth):
                                    CentralGrid(nCells, meshWidth)
{
};

void LbmDiscretization::calcRho(int i, int j, int k)
{
    double sum = 0.0;
        for(int l; l < 19 ; k++)
            {
                sum = sum + pdf(i,j,k,l);
            }
    rho(i,j, k) = sum;
};