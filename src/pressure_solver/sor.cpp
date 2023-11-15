#include "pressure_solver/sor.h"

SOR::SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega)
:PressureSolver(discretization, epsilon, maximumNumberOfIterations), omega_(omega)
{
}

void SOR::solve()
{
    double dxdx = pow(discretization_->dx(),2);
    double dydy = pow(discretization_->dy(),2);
    double factor_dx_dy =omega_* (dxdx*dydy)/(2*(dxdx+dydy)); 
    double eps2 = pow(epsilon_,2);
    const int  N = discretization_->nCells()[0]*discretization_->nCells()[1];

    double  residuum = 0.;
    int     iterations = 0;

    bool doSor = true;

    while(doSor)
    {
        for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++)
        {
            for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++)
            {
                    double dpdx = (discretization_->p(i-1,j) + discretization_->p(i+1,j))/dxdx;
                    double dpdy = (discretization_->p(i,j-1) + discretization_->p(i,j+1))/dydy; 

                    discretization_->p(i,j) = (1-omega_)*discretization_->p(i,j)+factor_dx_dy * (dpdx+dpdy - discretization_->rhs(i,j));
            } 
        }
        
        
        iterations++;
        setBoundaryValues();

        for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
            for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
                    double d2pdxx = (discretization_->p(i-1,j) - 2 * discretization_->p(i, j) + discretization_->p(i+1,j))/dxdx;
                    double d2pdyy = (discretization_->p(i,j-1) - 2 * discretization_->p(i, j) + discretization_->p(i,j+1))/dydy; 

                    residuum += pow(d2pdxx + d2pdyy - discretization_->rhs(i,j),2);
            } 
            residuum = residuum/N;
        }
            
        if(residuum < eps2)
        {
            doSor = false;
        }
        else if(iterations == maximumNumberOfIterations_)
        {
            doSor = false;
        }
    }
}