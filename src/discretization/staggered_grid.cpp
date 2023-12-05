#include "discretization/staggered_grid.h"
#include "iostream"

// construct staggered grid 
StaggeredGrid::StaggeredGrid(std::shared_ptr<Partitioning> partitioning, std::array< double, 2 > meshWidth):
    partitioning_(partitioning),
    nCells_(partitioning->getNCellsLocal()),    
    meshWidth_(meshWidth),
        f_({nCells_[0]+2, nCells_[1]+2}, {meshWidth_[0],        meshWidth_[1]/2.0}, meshWidth),
        g_({nCells_[0]+2, nCells_[1]+2}, {meshWidth_[0]/2.0,    meshWidth_[1]},     meshWidth),
        p_({nCells_[0]+2, nCells_[1]+2}, {meshWidth_[0]/2.0,    meshWidth_[1]/2.0}, meshWidth),
        u_({nCells_[0]+2, nCells_[1]+2}, {meshWidth_[0],        meshWidth_[1]/2.0}, meshWidth),
        v_({nCells_[0]+2, nCells_[1]+2}, {meshWidth_[0]/2.0,    meshWidth_[1]},     meshWidth),
        rhs_({nCells_[0]+2, nCells_[1]+2},           {meshWidth_[0]/2.0,    meshWidth_[1]/2.0}, meshWidth)

    {
};

//get the mesh width, i.e. the length of a single cell in x and y direction 
const std::array<double, 2> StaggeredGrid::meshWidth() const
{
    return meshWidth_;
}

//get number of cells in each coordinate direction 
const std::array<int, 2> StaggeredGrid::nCells() const
{
    return nCells_;
}

//get a reference to field variable u 
const FieldVariable & StaggeredGrid::p()const
{
    return p_;
}

//get the mesh width in x-direction, δx 
double StaggeredGrid::dx() const
{
    return meshWidth_[0];
}

//get the mesh width in y-direction, δy
double StaggeredGrid::dy() const
{
    return meshWidth_[1];
}

//access value of F in element (i,j) 
double & StaggeredGrid::f(int i, int j)
{
    return f_(i-uIBegin(),j-uJBegin());
}

//access value of G in element (i,j) 
double & StaggeredGrid::g(int i, int j)
{
    return g_(i-vIBegin(),j-vJBegin());
}

//access value of p in element (i,j) 
double StaggeredGrid::p(int i, int j) const
{
    return p_(i-pIBegin(),j-pJBegin());
}

//access value of p in element (i,j) 
double &StaggeredGrid::p(int i, int j)
{
    return p_(i-pIBegin(),j-pJBegin());
}

//first valid index for p in x direction 
int StaggeredGrid::pIBegin() const
{
    return -1;
}

//last valid index for p in x direction 
int StaggeredGrid::pIEnd() const
{
    return nCells_[0]+1;
}

//first valid index for p in y direction 
int StaggeredGrid::pJBegin() const
{
    return -1;
}

//last valid index for p in y direction 
int StaggeredGrid::pJEnd() const
{
    return nCells_[1]+1;
}

//access value of rhs in element (i,j) 
double & StaggeredGrid::rhs(int i, int j)
{
    return rhs_(i+1,j+1);
}

//get a reference to field variable u 
const FieldVariable & StaggeredGrid::u() const
{
    return u_;
}

//access value of u in element (i,j) 
double StaggeredGrid::u(int i, int j) const
{
    return u_(i-uIBegin(),j-uJBegin());
}

//access value of u in element (i,j) 
double & StaggeredGrid::u(int i, int j)
{
   return u_(i-uIBegin(),j-uJBegin());
}

//first valid index for u in x direction 
int StaggeredGrid::uIBegin() const
{
    if(partitioning_->ownPartitionContainsLeftBoundary())
        return -1;
    return -2;
}

//last valid index for u in x direction 
int StaggeredGrid::uIEnd() const
{
    if(partitioning_->ownPartitionContainsRightBoundary())
        return nCells_[0];
    return nCells_[0]+1;
}

//first valid index for u in y direction 
int StaggeredGrid::uJBegin() const
{
    return -1;
}

//last valid index for u in y direction 
int StaggeredGrid::uJEnd() const
{
    return nCells_[1]+1;
}

//get a reference to field variable v
const FieldVariable & StaggeredGrid::v() const
{
    return v_;
}

//access value of v in element (i,j) 
double StaggeredGrid::v(int i, int j) const
{
    return v_(i-vIBegin(),j-vJBegin());
}

//access value of v in element (i,j) 
double & StaggeredGrid::v(int i, int j) 
{
    return v_(i-vIBegin(),j-vJBegin());
}

//first valid index for v in x direction 
int StaggeredGrid::vIBegin() const
{
    return -1;
}

//last valid index for v in x direction 
int StaggeredGrid::vIEnd() const
{
    return nCells_[0]+1;
}

//first valid index for v in y direction 
int StaggeredGrid::vJBegin() const
{
    if(partitioning_->ownPartitionContainsBottomBoundary())
        return -1;
    return -2;
}

//last valid index for v in y direction 
int StaggeredGrid::vJEnd() const
{   
    if(partitioning_->ownPartitionContainsTopBoundary())
        return nCells_[1];
    return nCells_[1]+1;
}
