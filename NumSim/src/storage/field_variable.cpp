#include "field_variable.h"
#include "array2d.h"
#include <tgmath.h>
#include "settings/settings.h"

FieldVariable::FieldVariable(std::array<int,2> size, std::array<double, 2> origin, std::array<double, 2> meshWidth)
: Array2D(size), origin_(origin), meshWidth_(meshWidth)
{
}

double FieldVariable::interpolateAt(double x, double y) const {
    
    int i = floor((x-origin_[0])/meshWidth_[0]);
    int j = floor((y-origin_[1])/meshWidth_[1]);
    
    double left_value = (*this)(i,j);
    double right_value = (*this)(i+1,j);
    
    double left_pos = i*meshWidth_[0];
    double right_pos = j*meshWidth_[1];
    
    double interp_value = left_value + (right_value-left_value)/(right_pos - left_pos)*(x-left_pos);
    
    return interp_value;
}