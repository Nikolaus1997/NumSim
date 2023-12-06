#include "discretization/central_differences.h"
#include <cmath>

CentralDifferences::CentralDifferences(std::shared_ptr<Partitioning> partitioning, std::array<double,2> meshWidth)
: Discretization(partitioning, meshWidth)
{
};

//compute the 1st derivative ∂ u^2 / ∂x
double CentralDifferences::computeDu2Dx(int i, int j) const
{
    const double u_half       = (u(i+1,j)+u(i,j))/2.0;
    const double u_minus_half = (u(i,j)+u(i-1,j))/2.0;
    
    const double solution = (pow(u_half,2)-pow(u_minus_half,2))/dx();
    return solution;
}

//compute the 1st derivative ∂ v^2 / ∂x
double 	CentralDifferences::computeDv2Dy (int i, int j) const
{
    const double v_half       = (v(i,j+1)+v(i,j))/2.0;
    const double v_minus_half = (v(i,j)+v(i,j-1))/2.0;
    
    const double solution = (pow(v_half,2)-pow(v_minus_half,2))/dy();
    
    return solution;
}

//compute the 1st derivative ∂ (uv) / ∂y
double 	CentralDifferences::computeDuvDy (int i, int j) const
{
    const double u_half       = (u(i,j+1)+u(i,j))/2.0;
    const double u_minus_half = (u(i,j)+u(i,j-1))/2.0;

    const double v_half       = (v(i+1,j)+v(i,j))/2.0;
    const double v_minus_half = (v(i+1,j-1)+v(i,j-1))/2.0;
    
    const double solution = (v_half*u_half-v_minus_half*u_minus_half)/dy();
    
    return solution;
}

//compute the 1st derivative ∂ (uv) / ∂x
double CentralDifferences::computeDuvDx(int i, int j) const
{
    const double u_half       = (u(i,j)+u(i,j+1))/2.0;
    const double u_half_minus = (u(i-1,j)+u(i-1,j+1))/2.0;

    const double v_half       = (v(i,j)+v(i+1,j))/2.0;
    const double v_half_minus = (v(i,j)+v(i-1,j))/2.0;

    const double solution     = (u_half*v_half-u_half_minus*v_half_minus)/dx();
    
    return solution;
}