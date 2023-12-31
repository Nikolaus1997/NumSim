#pragma once

#include "storage/field_variable.h"
#include "partitioning/partitioning.h"

/**
 * The pressure p and velocities u and v are stored in a staggered grid.
 * This class provides functions to access p, u, v, f, g and rhs.
*/

class StaggeredGrid 
{
public:
    // construct staggered grid
    StaggeredGrid(std::shared_ptr<Partitioning> partitioning, std::array< double, 2 > meshWidth);

    //get the mesh width, i.e. the length of a single cell in x and y direction 
    const std::array<double, 2> meshWidth() const;

    //get number of cells in each coordinate direction 
    const std::array<int, 2> nCells() const;

    //get a reference to field variable u 
    const FieldVariable &u() const;
    //get a reference to field variable v
    const FieldVariable &v() const;
    //get a reference to field variable p 
    const FieldVariable &p() const;

    //access value of u in element (i,j) 
    double u(int i, int j) const;
    double &u(int i, int j);

    //access value of v in element (i,j) 
    double v(int i, int j) const;
    double &v(int i, int j);

    //access value of p in element (i,j) 
    double p(int i, int j) const;
    double &p(int i, int j);

    //access value of rhs in element (i,j) 
    double &rhs(int i, int j);
    //access value of F in element (i,j) 
    double &f(int i, int j);
    //access value of G in element (i,j) 
    double &g(int i, int j);

    //get the mesh width in x-direction, δx 
    double dx() const;
    //get the mesh width in y-direction, δy
    double dy() const;

    //first valid index for u in x direction 
    int uIBegin() const;
    //last valid index for u in x direction 
    int uIEnd() const;
    //first valid index for u in y direction 
    int uJBegin() const;
    //last valid index for u in y direction 
    int uJEnd() const;

    //first valid index for u in x direction 
    int rhsIBegin() const;
    //last valid index for u in x direction 
    int rhsIEnd() const;
    //first valid index for u in y direction 
    int rhsJBegin() const;
    //last valid index for u in y direction 
    int rhsJEnd() const;

    //first valid index for v in x direction 
    int vIBegin() const;
    //last valid index for v in x direction 
    int vIEnd() const;
    //first valid index for v in y direction 
    int vJBegin() const;
    //last valid index for v in y direction 
    int vJEnd() const;

    //first valid index for p in x direction 
    int pIBegin() const;
    //last valid index for p in x direction 
    int pIEnd() const;
    //first valid index for p in y direction 
    int pJBegin() const;
    //last valid index for p in y direction 
    int pJEnd() const;

    std::array<int, 2> uSize() const;
    std::array<int, 2> vSize() const;
    protected:
        const std::array<int, 2> nCells_;
        const std::array<double, 2> meshWidth_;
        FieldVariable u_;
        FieldVariable v_;
        FieldVariable p_;
        FieldVariable rhs_;
        FieldVariable f_;
        FieldVariable g_;
        std::shared_ptr<Partitioning> partitioning_;
};