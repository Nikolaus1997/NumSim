#pragma once
#include "discretization/staggered_grid.h"

/**
 * This class implements the derivative calculations in the form of finite differences.
*/

class Discretization : public StaggeredGrid
{

  public:

  // construct the object with given number of cells in x and y direction
  Discretization(std::shared_ptr<Partitioning> partitioning, std::array<double,2> meshWidth);

  // compute the 1st derivative ∂ u^2 / ∂x
  virtual double computeDu2Dx(int i, int j) const = 0;

  // compute the 1st derivative ∂ v^2 / ∂x
  virtual double computeDv2Dy(int i, int j) const = 0;

  // compute the 1st derivative ∂ (uv) / ∂x
  virtual double computeDuvDx(int i, int j) const = 0;

  // compute the 1st derivative ∂ (uv) / ∂y
  virtual double computeDuvDy(int i, int j) const = 0;

  // compute the 2nd derivative ∂^2 u / ∂x^2
  virtual double 	computeD2uDx2 (int i, int j) const;

  // compute the 2nd derivative ∂^2 u / ∂y^2  
  virtual double 	computeD2uDy2 (int i, int j) const;

  // compute the 2nd derivative ∂^2 v / ∂x^2
  virtual double 	computeD2vDx2 (int i, int j) const;
  
  // compute the 2nd derivative ∂^2 v / ∂y^2 
  virtual double 	computeD2vDy2 (int i, int j) const;

  // compute the 1st derivative ∂p / ∂x
  virtual double 	computeDpDx (int i, int j) const;

  // compute the 1st derivative ∂p / ∂y
  virtual double 	computeDpDy (int i, int j) const;
};