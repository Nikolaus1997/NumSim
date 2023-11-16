#include "pressure_solver/pressure_solver.h"




PressureSolver::PressureSolver(std::shared_ptr< Discretization >  discretization, double epsilon, int maximumNumberOfIterations)
:discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations)
{
    N_= discretization_->nCells()[0]*discretization_->nCells()[1];
}

void PressureSolver::solve()
{
}

void PressureSolver::setBoundaryValues()
{
       
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
    {
        //p bottom  boundary condition
        discretization_->p(i,discretization_->pJBegin()) = discretization_->p(i,discretization_->pJBegin()+1);

        //p top boundary condition
        discretization_->p(i,discretization_->pJEnd()-1)= discretization_->p(i,discretization_->pJEnd()-2);
    }
    
    for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
    {
        //p left boundary condition
        discretization_->p(discretization_->pIBegin(),j) = discretization_->p(discretization_->pIBegin()+1, j);

        //p right boundary condition
        discretization_->p(discretization_->pIEnd()-1, j)= discretization_->p(discretization_->pIEnd()-2, j);
    }
}
void PressureSolver::computeResiduum()
{
    double dxdx = pow(discretization_->dx(),2);
    double dydy = pow(discretization_->dy(),2);
    residuum_ = 0.;
            for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
            for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
                    double d2pdxx = (discretization_->p(i-1,j) - 2 * discretization_->p(i, j) + discretization_->p(i+1,j))/dxdx;
                    double d2pdyy = (discretization_->p(i,j-1) - 2 * discretization_->p(i, j) + discretization_->p(i,j+1))/dydy; 

                    residuum_ += pow(d2pdxx + d2pdyy - discretization_->rhs(i,j),2);
            } 
            
        }
}