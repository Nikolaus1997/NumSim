#include "computation.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/sor.h"

void Computation::initialize(int argc, char *argv[])
{
    //TODO
}

void Computation::runSimulation()
{
    //TODO
}

void Computation::computeTimeStepWidth()
{
    double dt_diffusion = settings_.re/2*(pow(discretization_->dx(),2)*pow(discretization_->dy(),2))/(pow(discretization_->dx(),2)+pow(discretization_->dy(),2));
    
    double dt_u = discretization_->dx()/discretization_->u().absMax();

    double dt_v = discretization_->dy()/discretization_->v().absMax();

    dt_ = settings_.tau*std::min(dt_diffusion,dt_u,dt_v);
}

void Computation::applyBoundaryValues()
{
    //TODO
}

void Computation::computePreliminaryVelocities()
{

    for (int i = discretization_->uIBegin()+1; i < discretization_->uIEnd(); i++)
    {
        for (int j = discretization_->uJBegin()+1; j < discretization_->uJEnd(); j++)
        {   
            double diffusion_u       = discretization_->computeD2uDx2 + discretization_->computeD2uDy2;
            double convection_u      = discretization_->computeDu2Dx  + discretization_->computeDuvDy;
            
            discretization_->f(i,j) = u(i,j) + dt_*(diffusion_u/settings_.re - convection_u + settings_.g[0]);
        }
        
    }

    for (int i = discretization_->vIBegin()+1; i < discretization_-> vIEnd(); i++)
    {
        for (int j = discretization_->vJBegin()+1; j < discretization_-> vJEnd(); j++)
        {
            double diffusion_v  = discretization_->computeD2vDx2 + discretization_->computeD2vDy2;
            double convection_v = discretization_->computeDv2Dy  + discretization_->computeDuvDx;
            
            discretization_->g(i,j) = v(i,j) + dt_*(diffusion_v/settings_.re - convection_v + settings_.g[1]);
        }
    }
    
}

void Computation::computeRightHandSide()
{
    for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd(); i++)
    {
        for (int j = discretization_->pJBegin()+1; j < discretization_-> pJEnd(); j++)
        {
            double dfdx = (discretization_->f(i,j) - discretization_-> f(i-1,j))/discretization_->dx();
            double dgdy = (discretization_->g(i,j) - discretization_-> g(i,j-1))/discretization_->dy();

            discretization_->rhs(i,j) = 1/dt_*(dfdx + dgdy);
        }
    }
}

void Computation::computePressure()
{
    pressureSolver_->solve();
}

void Computation::computeVelocities()
{
    //TODO
}