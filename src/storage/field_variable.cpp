#include "storage/field_variable.h"
#include "storage/array2d.h"
#include <tgmath.h>
#include "settings/settings.h"

FieldVariable::FieldVariable(std::array<int,2> size, std::array<double, 2> origin, std::array<double, 2> meshWidth)
: Array2D(size), origin_(origin), meshWidth_(meshWidth)
{
}

double FieldVariable::interpolateAt(double x, double y) const {
    
    int i = floor((x-origin_[0])/meshWidth_[0])+1;
    int j = floor((y-origin_[1])/meshWidth_[1])+1;
    
    if(i==size_[0]-1)
        i--;
    
    if(j==size_[1]-1)
        j--;

    double left         =       (*this)(i,j);
    double right        =       (*this)(i+1,j);
    double upper_left   =       (*this)(i,j+1);
    double upper_right  =       (*this)(i+1,j+1);

    double x_left   = i*meshWidth_[0] + origin_[0];
    double x_right   = meshWidth_[0] + x_left;

    double y_left   = j*meshWidth_[1] + origin_[1];
    double y_right   = meshWidth_[1] + y_left;
    
    double const interp_value = 1.0/(meshWidth_[0]*meshWidth_[1])*(
                                                                    left*(x_right-x)*(y_right-y)+
                                                                    right*(x-x_left)*(y_right-y)+
                                                                    upper_left*(x_right-x)*(y-y_left)+
                                                                    upper_right*(x-x_left)*(y-y_left)
        );
    
    return interp_value;
}
double FieldVariable::maxAbs()
{
    return 0.0;
}
