#include "storage/field_variable.h"
#include "staggered_grid.h"

StaggeredGrid::StaggeredGrid(std::array< int, 2 > nCells, std::array< double, 2 > meshWidth):
    nCells_(nCells),    
    meshWidth_(meshWidth), 
        FieldVariable u_    ({uIEnd()-uIBegin(), uJEnd()-uJBegin()}, {meshWidth_[0],        meshWidth_[1]/2.0}, meshWidth),
        FieldVariable v_    ({vIEnd()-vIBegin(), vJEnd()-vJBegin()}, {meshWidth_[0]/2.0,    meshWidth_[1]},     meshWidth),
        FieldVariable p_    ({pIEnd()-pIBegin(), pJEnd()-pJBegin()}, {meshWidth_[0]/2.0,    meshWidth_[1]/2.0}, meshWidth),
        FieldVariable rhs_  ({nCells_[0]+2, nCells_[1]+2},           {meshWidth_[0]/2.0,    meshWidth_[1]/2.0}, meshWidth),
        FieldVariable f_    ({uIEnd()-uIBegin(), uJEnd()-uJBegin()}, {meshWidth_[0],        meshWidth_[1]/2.0}, meshWidth),
        FieldVariable g_    ({vIEnd()-vIBegin(), vJEnd()-vJBegin()}, {meshWidth_[0]/2.0,    meshWidth_[1]},     meshWidth)
{
}

const std::array<double, 2> StaggeredGrid::meshWidth() const
{
    return meshWidth_;
}

const std::array<int, 2> StaggeredGrid::nCells() const
{
    return nCells_;
}

const FieldVariable & StaggeredGrid::p()const
{
    return p_;
}

double StaggeredGrid::dx() const
{
    return meshWidth_[0];
}

double StaggeredGrid::dy() const
{
    return meshWidth_[1];
}

double & StaggeredGrid::f(int i, int j)
{
    return f_(i,j);
}

double & StaggeredGrid::g(int i, int j)
{
    return g_(i,j);
}

double StaggeredGrid::p(int i, int j) const
{
    return p_(i,j);
}

double &StaggeredGrid::p(int i, int j)
{
    return p_(i,j);
}

int StaggeredGrid::pIBegin() const
{
    return -1;
}

int StaggeredGrid::pIEnd() const
{
    return nCells[0];
}

int StaggeredGrid::pJBegin() const
{
    return -1;
}

int StaggeredGrid::pJEnd() const
{
    return nCells[1];
}

double & StaggeredGrid::rhs(int i, int j)
{
    return rhs_(i,j);
}

const FieldVariable & StaggeredGrid::u() const
{
    return u_;
}

double StaggeredGrid::u(int i, int j) const
{
    return u_(i,j);
}

double StaggeredGrid::u(int i, int j)
{
   return u_(i,j);
}

int StaggeredGrid::uIBegin() const
{
    return -1;
}

int StaggeredGrid::uIEnd() const
{
    return nCells[0]-1;
}

int StaggeredGrid::uJBegin() const
{
    return -1;
}

int StaggeredGrid::uJEnd() const
{
    return nCells[1];
}

const FieldVariable & StaggeredGrid::v() const
{
    return v_;
}

double StaggeredGrid::v(int i, int j) const
{
    return v_(i,j);
}

double & StaggeredGrid::v(int i, int j) 
{
    return v_(i,j);
}

int StaggeredGrid::vIBegin() const
{
    return -1;
}

int StaggeredGrid::vIEnd() const
{
    return nCells[0];
}

int StaggeredGrid::vJBegin() const
{
    return -1;
}

int StaggeredGrid::vJEnd() const
{
    return nCells[1]-1;
}