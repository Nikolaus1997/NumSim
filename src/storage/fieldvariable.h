#pragma once

#include "storage/array3d.h"

/**
 * A field variable is the discretization of a scalar function f(x) with x in the computational domain.
 *  More specifically, a scalar value is stored at discrete nodes/points. The nodes are arranged in an equidistant mesh
 *  with specified mesh width.
 * @param size: number of cells
 * @param origin: origin of the coordinate system
 * @param meshWidth: width of cells in both directions
 */

class FieldVariable : public Array3D
{
public:
    /*
     * constructor
     */
    FieldVariable(std::array<int, 3> size,
                  std::array<double, 3> origin,
                  std::array<double, 3> meshWidth);
    
    // interpolate at arbitrary position
    double interpolateAt(double x, double y) const;

    double absMax() const;

private:
    std::array<double, 3> meshWidth_;
    std::array<double, 3> origin_;
};
