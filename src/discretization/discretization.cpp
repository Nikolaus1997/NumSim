#include "discretization/discretization.h"

//! construct the object with given number of cells in x and y direction
Discretization::Discretization(std::array<int, 2> nCells, std::array<double, 2> meshWidth):
    StaggeredGrid(nCells, meshWidth)
{
}

// compute the 2nd derivative ∂^2 u / ∂x^2
double Discretization::computeD2uDx2(int i, int j) const
{
    return (u(i+1,j)-2*u(i,j)+u(i-1,j))/(dx()*dx());
}

// compute the 1st derivative ∂ u^2 / ∂x
double Discretization::computeD2uDy2(int i, int j) const
{
    return (u(i,j+1)-2*u(i,j)+u(i,j-1))/(dy()*dy());
}

// compute the 1st derivative ∂ v^2 / ∂x
double Discretization::computeD2vDx2(int i, int j) const
{
    return (v(i+1,j)-2*v(i,j)+v(i-1,j))/(dx()*dx());
}

// compute the 2nd derivative ∂^2 v / ∂y^2 
double Discretization::computeD2vDy2(int i, int j) const
{
    return (v(i,j+1)-2*v(i,j)+v(i,j-1))/(dy()*dy());
}

// compute the 1st derivative ∂p / ∂x
double Discretization::computeDpDx(int i, int j) const
{
    return (p(i+1,j)-p(i,j))/dx();
}

// compute the 1st derivative ∂p / ∂y
double Discretization::computeDpDy(int i, int j) const
{
    return (p(i,j+1)-p(i,j))/dy();
}