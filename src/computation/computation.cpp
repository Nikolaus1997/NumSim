#include "computation.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/sor.h"
#include "settings/settings.h"

void Computation::initialize(std::string filename)
{
    settings_.loadFromFile(filename);
    settings_.printSettings();

    meshWidth_[0] = settings_.physicalSize[0]/settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1]/settings_.nCells[1];

    std::cout<<settings_.nCells[0]<<std::endl;
    
    std::cout<<"==========================================================================================================================="<<std::endl;
    if (settings_.useDonorCell) {
        std::cout<<"Using DonorCell..."<<std::endl;
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    }
    else {
        std::cout<<"Using CentralDifferences..."<<std::endl;
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    }
    
    if (settings_.pressureSolver == "SOR") {
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon,
                                                settings_.maximumNumberOfIterations, settings_.omega);
    }
    else if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon,
                                                settings_.maximumNumberOfIterations);

    } else {
        std::cout << "Solver not found!" << std::endl;
    }

    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
}

void Computation::runSimulation()
{
    int t_i = 0;
    double time = 0.0;
    while(time<settings_.endTime){
        t_i++;
        //std::cout << "Before applyBoundaryValues: " << discretization_->f(1,1) << std::endl;
        applyBoundaryValues();
        //std::cout << "After applyBoundaryValues: " << discretization_->f(1,1) << std::endl;
        applyBoundaryValuesFandG();

        computeTimeStepWidth();

        std::cout << "Time: " << time << ", dt: " << dt_ << std::endl;
        if (time + dt_>settings_.endTime){
            dt_ = settings_.endTime-time;
        }
        time = time + dt_;

        computePreliminaryVelocities();

        computeRightHandSide();

         computePressure();

         computeVelocities();

        //outputWriterParaview_->writeFile(time);
    }
}

void Computation::computeTimeStepWidth()
{
    double dt_diffusion = settings_.re/2*(pow(discretization_->dx(),2)*pow(discretization_->dy(),2))/(pow(discretization_->dx(),2)+pow(discretization_->dy(),2));
    
    double max_abs_u = 0.0;
    double max_abs_v = 0.0;
    double dt_u = 0.0;
    double dt_v = 0.0;
    double min_vel = 0.0;

    for(int i = discretization_->uIBegin(); i <=discretization_->uIEnd() ;i++){
        for(int j = discretization_->uJBegin(); j <= discretization_->uJEnd() ;j++){
            if(abs(discretization_->u(i,j)) > max_abs_u){
                max_abs_u = abs(discretization_->u(i,j));
            }
        }
    }

    for(int i = discretization_->vIBegin()+1; i <=discretization_->vIEnd() ;i++){
        for(int j = discretization_->vJBegin()+1; j <= discretization_->vJEnd() ;j++){
            if(abs(discretization_->v(i,j)) > max_abs_v){
                max_abs_v = abs(discretization_->v(i,j));
            }
    }
    }
    if(max_abs_u!=0.0){
        dt_u = discretization_->dx()/max_abs_u;
    }

    if(max_abs_v!=0.0){
        dt_v = discretization_->dy()/max_abs_v;
    }

    min_vel = std::min(dt_u,dt_v);
    if(min_vel == 0.0){
        dt_ = dt_diffusion*settings_.tau;
    }else{
        dt_ = settings_.tau*std::min(dt_diffusion,min_vel);
    }
}

void Computation::applyBoundaryValues()
{

    for(int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); i++)
    {
        //bottom boundary conditions u
        discretization_->u(i,discretization_->uJBegin()) = 2*settings_.dirichletBcBottom[0]-discretization_->u(i,discretization_->uJBegin()+1);
    

        //top boundary conditions u
        discretization_->u(i,discretization_->uJEnd()) = 2*settings_.dirichletBcTop[0]-discretization_->u(i,discretization_->uJEnd()-1);

    }

    for(int j = discretization_->uJBegin(); j <=discretization_->uJEnd();j++)
    {
        //left boundary conditions u
        discretization_->u(discretization_->uIBegin(),j) = settings_.dirichletBcLeft[0];

        //right boundary conditions u
        discretization_->u(discretization_->uIEnd(),j) = settings_.dirichletBcRight[0];
    }

    for(int i = discretization_->vIBegin(); i <= discretization_->vIEnd();i++)
    {
        //top boundary conditions v
        discretization_->v(i,discretization_->vJEnd()) = settings_.dirichletBcTop[1];

        //bottom boundary conditions v
        discretization_->v(i, discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
    }

    for(int j = discretization_->vJBegin(); j <= discretization_->vJEnd(); j++)
    {
        //left boundary condition
        discretization_->v(discretization_->vIBegin(),j) = 2.0*settings_.dirichletBcLeft[1]-discretization_->v(discretization_->vIBegin()+1,j);

        //right boundary condition
        discretization_->v(discretization_->vIEnd(),j)  = 2.0*settings_.dirichletBcRight[1]-discretization_->v(discretization_->vIEnd()-1,j);
    }

}

void Computation::applyBoundaryValuesFandG(){
    for (int i = discretization_->uIBegin(); i <=discretization_->uIEnd(); i++)
    {
        //top boundary conditions for u
        discretization_->f(i, discretization_->uJEnd()) = discretization_->u(i,discretization_->uJEnd());

        //bottom boundary conditions for u
        discretization_->f(i, discretization_->uJBegin()) = discretization_->u(i,discretization_->uJBegin());
    }

    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++)
    {   
        //left boundary condition for u
        discretization_->f(discretization_->uIBegin(), j) = discretization_->u(discretization_->uIBegin(),j);

        //right boundary condition for u
        discretization_->f(discretization_->uIEnd(), j) = discretization_->u(discretization_->uIEnd(), j);
    }

    for (int i = discretization_->vIBegin(); i <=discretization_->vIEnd(); i++)
    {
        //top boundary conditions for v
        discretization_->g(i, discretization_->vJEnd()) = discretization_->v(i,discretization_->vJEnd());

        //bottom boundary conditions for v
        discretization_->g(i, discretization_->vJBegin()) = discretization_->v(i,discretization_->vJBegin());
    }
    
    for (int j = discretization_->vJBegin(); j <= discretization_->vJEnd(); j++)
    {
        //left boundary conditions for v
        discretization_->g(discretization_->vIBegin(),j) = discretization_->v(discretization_->vIBegin(),j);

        //right boundary condtions for v
        discretization_->g(discretization_->vIEnd(),j) = discretization_->v(discretization_->vIEnd(),j);
    }

}

void Computation::computePreliminaryVelocities()
{

    for (int i = discretization_->uIBegin()+1; i < discretization_->uIEnd(); i++)
    {
        for (int j = discretization_->uJBegin()+1; j < discretization_->uJEnd(); j++)
        {   
            double diffusion_u       = discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j);
            double convection_u      = discretization_->computeDu2Dx(i,j)  + discretization_->computeDuvDy(i,j);
            
            discretization_->f(i,j) = discretization_->u(i,j) + dt_*(diffusion_u/settings_.re - convection_u + settings_.g[0]);
        }
        
    }

    for (int i = discretization_->vIBegin()+1; i < discretization_-> vIEnd(); i++)
    {
        for (int j = discretization_->vJBegin()+1; j < discretization_-> vJEnd(); j++)
        {
            double diffusion_v  = discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j);
            double convection_v = discretization_->computeDv2Dy(i,j)  + discretization_->computeDuvDx(i,j);
            
            discretization_->g(i,j) = discretization_->v(i,j) + dt_*(diffusion_v/settings_.re - convection_v + settings_.g[1]);
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
    //set values for u
    for(int i = discretization_->uIBegin()+1; i < discretization_->uIEnd(); i++){
        for(int j = discretization_->uJBegin()+1; j < discretization_->uJEnd(); j++){
            discretization_->u(i,j) = discretization_->f(i,j)-dt_*discretization_->computeDpDx(i,j);
        }
    }

    //set values for v
    for(int i = discretization_->vIBegin()+1; i < discretization_->vIEnd(); i++){
        for(int j = discretization_->vJBegin()+1; j < discretization_-> vJEnd(); j++){
            discretization_->v(i,j) = discretization_->g(i,j)-dt_*discretization_->computeDpDy(i,j);
        }
    }
}

