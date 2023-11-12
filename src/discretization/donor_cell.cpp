#include "donor_cell.h"
#include <cmath>

DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth, double alpha): 
Discretization(nCells, meshWidth), alpha_(alpha)
{

};

double DonorCell::computeDu2Dx(int i, int j) const
{
    const double u_first  = (u(i,j)+u(i+1,j))/2.0;
    const double u_second = (u(i-1,j)+u(i,j))/2.0;

    const double u_first_minus  = (u(i,j)-u(i+1,j))/2.0;
    const double u_second_minus = (u(i-1,j)-u(i,j))/2.0;
    
    const double first_term = pow(u_first,2) - pow(u_second,2);
    const double secod_term = abs(u_first)*u_first_minus - abs(u_second)*u_second_minus;

    const double solution   =  1/dx() * (first_term + alpha_*secod_term);

    return solution;
}

double DonorCell::computeDv2Dy(int i, int j) const
{
    const double v_first  = (v(i,j)+v(i,j+1))/2.0;
    const double v_second = (v(i,j-1)+v(i,j))/2.0;

    const double v_first_minus  = (v(i,j)-v(i,j+1))/2.0;
    const double v_second_minus = (v(i,j-1)-v(i,j))/2.0;
    
    const double first_term = pow(v_first,2) - pow(v_second,2);
    const double secod_term = abs(v_first)*v_first_minus - abs(v_second)*v_second_minus;

    const double solution   =  1/dy() * (first_term + alpha_*secod_term);

    return solution;
}

double DonorCell::computeDuvDy(int i, int j) const
{
    const double v_first  = (v(i,j)+v(i+1,j))/2.0;
    const double u_second = (u(i,j)+u(i,j+1))/2.0;

    const double v_first_minus  = (v(i,j-1)+v(i+1,j-1))/2.0;
    const double u_second_minus = (u(i,j-1)+u(i,j))/2.0;

    const double u_first_minus_minus  = (u(i,j)-u(i,j+1))/2.0;
    const double u_second_minus_minus = (u(i,j-1)-u(i,j))/2.0; 
    
    const double first_term = v_first*u_second - v_first_minus*u_second_minus;
    const double secod_term = abs(v_first)*u_first_minus_minus - abs(v_first_minus)*u_second_minus_minus;

    const double solution   =  1/dy() * (first_term + alpha_*secod_term);

    return solution;
}

double DonorCell::computeDuvDx(int i, int j) const
{
    const double u_first    = (u(i,j)+u(i,j+1))/2.0;
    const double v_second   = (v(i,j)+v(i+1,j))/2.0;

    const double u_second_minus = (u(i-1,j)+u(i-1,j+1))/2.0;
    const double v_second_minus = (v(i-1,j)+v(i,j))/2.0;

    const double v_second_minus_minus       = (v(i,j)-v(i+1,j))/2.0;
    const double v_second_minus_minus_minus = (v(i-1,j)-v(i,j))/2.0;

    const double first_term  = u_first*v_second-u_second_minus*v_second_minus;
    const double second_term = abs(u_first)*v_second_minus_minus-abs(u_second_minus*v_second_minus_minus_minus);
    
    const double solution = 1/dx()*(first_term + alpha_*second_term);

    return solution;
}