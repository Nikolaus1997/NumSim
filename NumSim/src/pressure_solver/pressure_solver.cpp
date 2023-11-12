#include "pressure_solver.h"



PressureSolver::PressureSolver(std::shared_ptr< Discretization >  discretization, double epsilon, int maximumNumberOfIterations)
:discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations)
{

}

void PressureSolver::setBoundaryValues()
{
    
    for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++)
    {
        //p left boundary condition
        discretization_->p(discretization_->pIBegin(),j) = discretization_->p(discretization_->pIBegin()+1, j);

        //p right boundary condition
        discretization_->p(discretization_->pIEnd(), j)= discretization_->p(discretization_->pIEnd()-1, j);
    }
    
    for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++)
    {
        //p bottom  boundary condition
        discretization_->p(i,discretization_->pJBegin()) = discretization_->p(i,discretization_->pJBegin()+1);

        //p top boundary condition
        discretization_->p(i,discretization_->pJEnd())= discretization_->p(i,discretization_->pJEnd()-1);
    }
    
}