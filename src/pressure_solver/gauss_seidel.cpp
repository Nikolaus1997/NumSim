#include "gauss_seidel.h"

GaussSeidel::GaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations)
: PressureSolver(discretization, epsilon, maximumNumberOfIterations)
{
}

void GaussSeidel::solve()
{
    double dxdx = pow(discretization_->dx(),2);
    double dydy = pow(discretization_->dy(),2);
    double factor_dx_dy = dxdx*dydy/(2(dxdx+dydy));    

    double  residuum = 0.;
    int     iterations = 0;

    bool doGaussSeidel = true;

    while(doGaussSeidel == true)
    {
        for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd(); i++)
        {
            for (int j = discretization_->pJBegin()+1; j < discretization_->pJBegin(); j++)
            {
                    double dpdx = (discretization_->p(i-1,j) + discretization_->p(i+1,j))/dxdx;
                    double dpdy = (discretization_->p(i,j-1) + discretization_->p(i,j+1))/dydy; 

                    discretization_->p(i,j) = factor_dx_dy * (dpdx+dpdy - discretization_->rhs(i,j));
            } 
        }
        
        
        iterations++;
        setBoundaryValues();

        for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd(); i++)
        {
            for (int j = discretization_->pJBegin()+1; j < discretization_->pJBegin(); j++)
            {
                    double d2pdxx = (discretization_->p(i-1,j) - 2 * discretization_->p(i, j) + discretization_->p(i+1,j))/dxdx;
                    double d2pdyy = (discretization_->p(i,j-1) - 2 * discretization_->p(i, j) + discretization_->p(i,j+1))/dydy; 

                    residuum = std::max(std::abs(d2pdxx+d2pdyy - discretization_->rhs(i,j)), residuum);
            } 
        }
            
        if(residuum < epsilon_)
        {
            doGaussSeidel = false;
        }
        else if(iterations = maximumNumberOfIterations_)
        {
            doGaussSeidel = false;
        }
    }
}
