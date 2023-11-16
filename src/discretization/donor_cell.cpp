#include "discretization/donor_cell.h"
#include <cmath>

//constructor
DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth, double alpha): 
Discretization(nCells, meshWidth), alpha_(alpha)
{
};

//compute the 1st derivative ∂ u^2 / ∂x
double DonorCell::computeDu2Dx(int i, int j) const
{
    const double u_first  = (u(i,j)+u(i+1,j))/2.0;
    const double u_second = (u(i-1,j)+u(i,j))/2.0;

    const double u_first_minus  = (u(i,j)-u(i+1,j))/2.0;
    const double u_second_minus = (u(i-1,j)-u(i,j))/2.0;
    
    const double first_term = (pow(u_first,2) - pow(u_second,2))/dx();
    const double secod_term = (fabs(u_first)*u_first_minus - fabs(u_second)*u_second_minus)/dx();

    const double solution   =  first_term + alpha_*secod_term;

    return solution;

}

//compute the 1st derivative ∂ v^2 / ∂x 
double DonorCell::computeDv2Dy(int i, int j) const
{
    const double v_first  = (v(i,j)+v(i,j+1))/2.0;
    const double v_second = (v(i,j-1)+v(i,j))/2.0;

    const double v_first_minus  = (v(i,j)-v(i,j+1))/2.0;
    const double v_second_minus = (v(i,j-1)-v(i,j))/2.0;
    
    const double first_term = (pow(v_first,2) - pow(v_second,2))/dy();
    const double secod_term = (fabs(v_first)*v_first_minus - fabs(v_second)*v_second_minus)/dy();

    const double solution   =  first_term + alpha_*secod_term;

    return solution;

}

//compute the 1st derivative ∂ (uv) / ∂y 
double DonorCell::computeDuvDy(int i, int j) const
{
    const double v_first  = (v(i,j)+v(i+1,j))/2.0;
    const double u_second = (u(i,j)+u(i,j+1))/2.0;

    const double v_first_minus  = (v(i,j-1)+v(i+1,j-1))/2.0;
    const double u_second_minus = (u(i,j-1)+u(i,j))/2.0;

    const double u_first_minus_minus  = (u(i,j)-u(i,j+1))/2.0;
    const double u_second_minus_minus = (u(i,j-1)-u(i,j))/2.0; 
    
    const double first_term = (v_first*u_second - v_first_minus*u_second_minus) /dy();
    const double secod_term = (fabs(v_first)*u_first_minus_minus - fabs(v_first_minus)*u_second_minus_minus)/dy();

    const double solution   =  first_term + alpha_*secod_term;

    return solution;
}

//compute the 1st derivative ∂ (uv) / ∂x
double DonorCell::computeDuvDx(int i, int j) const {
    const double u_donor_right= (u(i,j+1) + u(i,j)) / 2.0;
    const double u_donor_left =  (u(i-1,j+1) + u(i-1,j)) / 2.0;

    const double v_donor_right = (v(i,j) + v(i+1,j)) / 2.0;
    const double v_donor_left = (v(i-1,j) + v(i,j)) / 2.0;

    const double v_minus_donor_right = (v(i,j) - v(i+1,j)) / 2.0;
    const double v_minus_donor_left = (v(i-1,j) - v(i,j)) / 2.0;

    const double first_term = (u_donor_right* v_donor_right - u_donor_left * v_donor_left) / dx();
    const double second_term = (fabs(u_donor_right) * v_minus_donor_right - fabs(u_donor_left) * v_minus_donor_left) / dx();

    const double solution = first_term + alpha_ * second_term;
    
    return solution;
}