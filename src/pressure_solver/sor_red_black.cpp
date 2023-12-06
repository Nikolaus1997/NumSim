#include "pressure_solver/sor_red_black.h"
#include "pressure_solver/pressure_solver_parallel.h"

SORRedBlack::SORRedBlack(std::shared_ptr<Discretization> discretization, 
                    double epsilon, 
                    int maximumNumberOfIterations, 
                    double omega,
                    std::shared_ptr<Partitioning> partitioning) :
    PressureSolverParallel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations),
    omega_(omega) {
}

void SORRedBlack::solve() {
    bool doSor = true;
    double dxdx = pow(discretization_->dx(),2);
    double dydy = pow(discretization_->dy(),2);
    double factor_omega = 1-omega_;
    double factor_dx_dy =omega_* (dxdx*dydy)/(2*(dxdx+dydy)); 
    double eps2 = pow(epsilon_,2);
    double  residuum_norm = 0.;
    int     iterations = 0;

    while(doSor){
        //black
        for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
            iStart = discretization_->pIBegin()+1
            for (int i = iStart; i < discretization_->pIEnd()-1; i+=2){
                    double dpdx = (discretization_->p(i-1,j) + discretization_->p(i+1,j))/dxdx;
                    double dpdy = (discretization_->p(i,j-1) + discretization_->p(i,j+1))/dydy; 
                    discretization_->p(i,j) = factor_omega*discretization_->p(i,j)+factor_dx_dy * (dpdx+dpdy - discretization_->rhs(i,j));
            } 
        }

        //communication

        //red
        for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
            iStart = discretization_->pIBegin()+2
            for (int i = iStart; i < discretization_->pIEnd()-1; i+=2){
                    double dpdx = (discretization_->p(i-1,j) + discretization_->p(i+1,j))/dxdx;
                    double dpdy = (discretization_->p(i,j-1) + discretization_->p(i,j+1))/dydy; 
                    discretization_->p(i,j) = factor_omega*discretization_->p(i,j)+factor_dx_dy * (dpdx+dpdy - discretization_->rhs(i,j));
            } 
        }

        //communication

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