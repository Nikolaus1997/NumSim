#include "storage/field_variable.h"
#include "storage/array2d.h"
#include <tgmath.h>
#include "settings/settings.h"
#include <assert.h>

FieldVariable::FieldVariable(std::array<int,2> size, std::array<double, 2> origin, std::array<double, 2> meshWidth)
: Array2D(size), origin_(origin), meshWidth_(meshWidth)
{
}

double FieldVariable::interpolateAt(double x, double y) const {
    //compute indices i and j of the cell at position x and y
    int i = floor((x-origin_[0])/meshWidth_[0])+1;
    int j = floor((y-origin_[1])/meshWidth_[1])+1;
    
    //if cell is at top or right boundary, move inside
    if(i==size_[0]-1){
        i--;
    }
    if(j==size_[1]-1){
        j--;
    }
    
    double left         =       (*this)(i,j);
    double right        =       (*this)(i+1,j);
    double upper_left   =       (*this)(i,j+1);
    double upper_right  =       (*this)(i+1,j+1);

    double x_left   = (i-1)*meshWidth_[0] + origin_[0];
    double x_right   = meshWidth_[0] + x_left;

    double y_left   = (j-1)*meshWidth_[1] + origin_[1];
    double y_right   = meshWidth_[1] + y_left;
    
    //bilinear interpolation formula
    double const interp_value = 1.0/(meshWidth_[0]*meshWidth_[1])*(
                                                                    left*(x_right-x)*(y_right-y)+
                                                                    right*(x-x_left)*(y_right-y)+
                                                                    upper_left*(x_right-x)*(y-y_left)+
                                                                    upper_right*(x-x_left)*(y-y_left)
        );
    
    return interp_value;

}
//TODO: implement this function in our way
double FieldVariable::interpolateAtParallel(double x, double y, std::shared_ptr<Partitioning> partitioning) const {
 {
    // Assert that the specified point is part of the domain
    #ifndef NDEBUG
    assert((0.0 <= x) && (x <= size_[0] * meshWidth_[0]));
    assert((0.0 <= y) && (y <= size_[1] * meshWidth_[1]));
    #endif

    // Determine i and j indices of the corresponding cell (shifted by origin)
    int i = (x - origin_[0]) / meshWidth_[0] + 1;
    int j = (y - origin_[1]) / meshWidth_[1] + 1;

    bool u = ((origin_[0] == meshWidth_[0]) && (origin_[1] == meshWidth_[1] / 2.));
    bool v = ((origin_[0] == meshWidth_[0] / 2.) && (origin_[1] == meshWidth_[1]));

    int i_test = i;
    int j_test = j;
    int x_size = size_[0];
    int y_size = size_[1];

    // Determines if GhostLayers exists and shifts the index if necessary and adjusts the boundary index
    if (!partitioning->ownPartitionContainsLeftBoundary() && u){
        i_test++;
        x_size--;
    }
    if (!partitioning->ownPartitionContainsRightBoundary() && u){
        x_size--;
    }
    if (!partitioning->ownPartitionContainsBottomBoundary() && v){
        j_test++;
        y_size--;
    }
    if (!partitioning->ownPartitionContainsTopBoundary() && v){
        y_size--;
    }

    // Special case: If we are on the upper of right boundary, use the cell in the interior
    if (i == x_size - 1)
        i--;
    if (j == y_size - 1)
        j--;
    // Special case: If we are on the upper of right boundary, use the cell in the interior
    if (i_test == size_[0] - 1)
        i_test--;
    if (j_test == size_[1] - 1)
        j_test--;


    // Obtain the values of the four neighbouring interpolation points
    double valueLeftBottom = (*this)(i_test, j_test);
    double valueLeftTop = (*this)(i_test, j_test + 1);
    double valueRightBottom = (*this)(i_test + 1, j_test);
    double valueRightTop = (*this)(i_test + 1, j_test + 1);

    // Determine the coordinates of the four neighbouring interpolation points
    double xLeft = origin_[0] + (i - 1) * meshWidth_[0];
    double xRight = xLeft + meshWidth_[0];
    double yBottom = origin_[1] + (j - 1) * meshWidth_[1];
    double yTop = yBottom + meshWidth_[1];

    /*
     * Bilinear interpolation:
     * We implement it as a repeated linear interpolation (first along x axis, then along y axis)
     * See also https://en.wikipedia.org/wiki/Bilinear_interpolation#Repeated_linear_interpolation
     */
    // 1) Use linear interpolation in x between left and right edge (each on the bottom and top edge)
    double interpBottom = ((xRight - x) * valueLeftBottom + (x - xLeft) * valueRightBottom) / (xRight - xLeft);
    double interpTop = ((xRight - x) * valueLeftTop + (x - xLeft) * valueRightTop) / (xRight - xLeft);
    // 2) Use linear interpolation in y between bottom and top edge
    double interp = ((yTop - y) * interpBottom + (y - yBottom) * interpTop) / (yTop - yBottom);

    return interp;
}


}
