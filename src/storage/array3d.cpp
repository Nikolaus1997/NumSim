#include "storage/array3d.h"
#include <iostream>
#include <cassert>

/**
 * This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 * @param size: number of cells
 */

Array3D::Array3D(std::array<int, 3> size):
        size_(size) {
    // allocate data, initialize to 0
    data_.resize(size_[0] * size_[1] * size_[2], 0.0);
}

/**
 * get number of cells
 * @return number of cells
 */
std::array<int, 3> Array3D::size() const {
    return size_;
}
/**
 * access the value at coordinate (i,j), declared not const, i.e. the value can be changed
 * @param i: discretized position in x direction
 * @param j: discretized position in y direction
 * @param k: discretized position in z direction
 * @return reference to value at the grid cell (i,j)
 */
double &Array3D::operator()(int i, int j, int k) {
    const int index = j * size_[0] + i + k*size_[0]*size_[1];

    // assert that indices are in range
    #ifndef NDEBUG
    assert((0 <= i) && (i < size_[0]));
    assert((0 <= j) && (j < size_[1]));
    assert((0 <= k) && (k < size_[2]));
    assert(j * size_[0] + i*size_[2] + k < (int) data_.size());
    #endif

    return data_[index];
}
/**
 * get the value at coordinate (i,j), declared const, i.e. it is not possible to change the value
 * @param i: discretized position in x direction
 * @param j: discretized position in y direction
 * @return value at the grid cell (i,j)
 */
double Array3D::operator()(int i, int j, int k) const {
    const int index = j * size_[0] + i + k*size_[0]*size_[1];

    // assert that indices are in range
    #ifndef NDEBUG
    assert((0 <= i) && (i < size_[0]));
    assert((0 <= j) && (j < size_[1]));
    assert((0 <= k) && (k < size_[2]));
    assert(j * size_[0] + i*size_[2] + k < (int) data_.size());
    #endif

    return data_[index];
}

/**
 * print out grid in two dimensions
 */

void Array3D::print() const {
    std::cout << std::endl << "----------" << std::endl;
    for (int j = size_[1] - 1; j >= 0; j--) {
        for (int i = 0; i < size_[0]; i++) {
            for (int k = 0; k < size_[1]; k++){
                std::cout << (*this)(i, j, k) << " | ";
            }
        }
        std::cout << std::endl << "----------" << std::endl;
    }
}