#include "pressure_solver/sor.h"

SOR::SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega)
:PressureSolver(discretization, epsilon, maximumNumberOfIterations), omega_(omega)
{
}

//solve the system of the Poisson equation for pressure using SOR
void SOR::solve()
{
    double dxdx = pow(discretization_->dx(),2);
    double dydy = pow(discretization_->dy(),2);
    double factor_omega = 1-omega_;
    double factor_dx_dy =omega_* (dxdx*dydy)/(2*(dxdx+dydy)); 
    double eps2 = pow(epsilon_,2);
    double  residuum_norm = 0.;
    int     iterations = 0;

    bool doSor = true;

    while(doSor){
        for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
            for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
                    double dpdx = (discretization_->p(i-1,j) + discretization_->p(i+1,j))/dxdx;
                    double dpdy = (discretization_->p(i,j-1) + discretization_->p(i,j+1))/dydy; 

                    discretization_->p(i,j) = factor_omega*discretization_->p(i,j)+factor_dx_dy * (dpdx+dpdy - discretization_->rhs(i,j));
            } 
        }
        
        iterations++;
        //set boundaries after each iteration
        setBoundaryValues();
        //compute residuum and check for convergence
        computeResiduum();  
        residuum_norm = residuum_/N_;
        if(residuum_norm < eps2 || iterations == maximumNumberOfIterations_)
        {
            doSor = false;
        }
    }
}