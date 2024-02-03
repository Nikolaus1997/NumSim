#include "storage/array4d.h"
#include <iostream>
#include <cassert>

/**
 * This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 * @param size: number of cells
 */

Array4D::Array4D(std::array<int, 4> size):
            size_(size) {
    // allocate data, initialize to 0
    data_.resize(size_[0] * size_[1] * size_[2]*size_[3], 0.0);
}

/**
 * get number of cells
 * @return number of cells
 */
std::array<int, 4> Array4D::size() const {
    return size_;
}
/**
 * access the value at coordinate (i,j), declared not const, i.e. the value can be changed
 * @param i: discretized position in x direction
 * @param j: discretized position in y direction
 * @param k: discretized position in z direction
 * @return reference to value at the grid cell (i,j)
 */
double &Array4D::operator()(int i, int j, int k, int l) {
    const int index = j * size_[0] + i + k*size_[0]*size_[1]+l*size_[0]*size_[1]*size_[2];

    // assert that indices are in range
    #ifndef NDEBUG
    assert((0 <= i) && (i < size_[0]));
    assert((0 <= j) && (j < size_[1]));
    assert((0 <= k) && (k < size_[2]));
    assert((0 <= l) && (l < size_[3]));    
    assert(j * size_[0] + i + k*size_[0]*size_[1]+l*size_[0]*size_[1]*size_[2]< (int) data_.size());
    #endif

    return data_[index];
}
/**
 * get the value at coordinate (i,j), declared const, i.e. it is not possible to change the value
 * @param i: discretized position in x direction
 * @param j: discretized position in y direction
 * @return value at the grid cell (i,j)
 */
double Array4D::operator()(int i, int j, int k, int l) const {
    const int index = j * size_[0] + i + k*size_[0]*size_[1]+l*size_[0]*size_[1]*size_[2];

    // assert that indices are in range
    #ifndef NDEBUG
    assert((0 <= i) && (i < size_[0]));
    assert((0 <= j) && (j < size_[1]));
    assert((0 <= k) && (k < size_[2]));
    assert((0 <= l) && (l < size_[3]));      
    assert(j * size_[0] + i + k*size_[0]*size_[1]+l*size_[0]*size_[1]*size_[2]< (int) data_.size());
    #endif

    return data_[index];
}

/**
 * print out grid in two dimensions
 */

void Array4D::print() const {
    std::cout << std::endl << "----------" << std::endl;
    for (int j = size_[1] - 1; j >= 0; j--) {
        for (int i = 0; i < size_[0]; i++) {
            for (int k = 0; k < size_[2]; k++){
                for (int l = 0; l < size_[3]; k++){
                std::cout << (*this)(i, j, k, l) << " | ";
            }
            }
        }
        std::cout << std::endl << "----------" << std::endl;
    }
}