#include "pressure_solver/gauss_seidel.h"

GaussSeidel::GaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations)
: PressureSolver(discretization, epsilon, maximumNumberOfIterations)
{
}

//solve the Poisson problem for the pressure, using the Gauss-Seidel method
void GaussSeidel::solve()
{
    double dxdx = pow(discretization_->dx(),2);
    double dydy = pow(discretization_->dy(),2);
    double factor_dx_dy = dxdx*dydy/(2*(dxdx+dydy));    
    double eps2 = epsilon_*epsilon_;
    double  residuum_norm = 0.;
    double  residuum = 0.;
    int     iterations = 0;

    bool doGaussSeidel = true;

    //Gauss Seidel method
    while(doGaussSeidel == true)
    {
        for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
            for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
                    double dpdx = (discretization_->p(i-1,j) + discretization_->p(i+1,j))/dxdx;
                    double dpdy = (discretization_->p(i,j-1) + discretization_->p(i,j+1))/dydy; 

                    discretization_->p(i,j) = factor_dx_dy * (dpdx+dpdy - discretization_->rhs(i,j));
            } 
        }
        iterations++;
        // set boundary values after each iteration
        setBoundaryValues();

        //calculate residuum and check for convergence
        computeResiduum();  
        residuum_norm = residuum_/N_;
        if(residuum_norm < eps2||iterations == maximumNumberOfIterations_)
        {
            doGaussSeidel = false;
        }
    }
}
