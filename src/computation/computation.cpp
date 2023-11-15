#include "computation.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/sor.h"
#include "settings/settings.h"

void Computation::initialize(std::string filename)
{   
    settings_ = Settings();
    settings_.loadFromFile(filename);
    settings_.printSettings();

    meshWidth_[0] = settings_.physicalSize[0]/settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1]/settings_.nCells[1];
    
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
    outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
}

void Computation::runSimulation()
{
    int t_i = 0;
    double time = 0.0;
    while(time<settings_.endTime){
        t_i++;
        applyBoundaryValues();
        applyBoundaryValuesFandG();;
        computeTimeStepWidth();

        if(dt_==0.0){
            std::cout<<"dt is zero, exiting..."<<std::endl; 
            break;
        }

        if (time + dt_>settings_.endTime){
            dt_ = settings_.endTime-time;
        }
        time = time + dt_;
        computePreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        outputWriterText_->writeFile(time);
        outputWriterParaview_->writeFile(time);
    }
}

void Computation::computeTimeStepWidth()
{
    double dt_diffusion = settings_.re / 2 / (1 / (discretization_->dx() * discretization_->dx()) + 1 / (discretization_->dy() * discretization_->dy()) );
    
    double max_abs_u = 0.0;
    double max_abs_v = 0.0;
    double dt_u = 0.0;
    double dt_v = 0.0;
    double vel_u = 0.0;
    double vel_v = 0.0;

    for(int i = discretization_->uIBegin(); i <discretization_->uIEnd() ;i++){
        for(int j = discretization_->uJBegin(); j < discretization_->uJEnd() ;j++){
            if(fabs(discretization_->u(i,j)) > max_abs_u){
                max_abs_u = fabs(discretization_->u(i,j));
            }
        }
    }

    for(int i = discretization_->vIBegin()+1; i <discretization_->vIEnd() ;i++){
        for(int j = discretization_->vJBegin()+1; j < discretization_->vJEnd() ;j++){
            if(fabs(discretization_->v(i,j)) > max_abs_v){
                max_abs_v = fabs(discretization_->v(i,j));
            }
    }
    }
    if(max_abs_u!=0){        
        vel_u = discretization_->dx()/max_abs_u;
        dt_= settings_.tau*std::min({vel_u, dt_diffusion});   
    }else if(max_abs_v!=0){
        vel_v = discretization_->dy()/max_abs_v;
        dt_= settings_.tau*std::min({vel_v, dt_diffusion}); 
    }else{    dt_= settings_.tau*std::min({dt_diffusion});
    }   


}

void Computation::applyBoundaryValues()
{

    for(int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        //bottom boundary conditions u
        discretization_->u(i,discretization_->uJBegin()) = 2.0*settings_.dirichletBcBottom[0]-discretization_->u(i,discretization_->uJBegin()+1);
    

        //top boundary conditions u
        discretization_->u(i,discretization_->uJEnd()-1) = 2.0*settings_.dirichletBcTop[0]-discretization_->u(i,discretization_->uJEnd()-2);

    }
    
    for(int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
    {
        //bottom boundary conditions v
        discretization_->v(i, discretization_->vJBegin()) = settings_.dirichletBcBottom[1];

        //top boundary conditions v
        discretization_->v(i,discretization_->vJEnd()-1) = settings_.dirichletBcTop[1];

    }

    for(int j = discretization_->uJBegin(); j <discretization_->uJEnd(); j++)
    {
        //left boundary conditions u
        discretization_->u(discretization_->uIBegin(),j) = settings_.dirichletBcLeft[0];

        //right boundary conditions u
        discretization_->u(discretization_->uIEnd()-1,j) = settings_.dirichletBcRight[0];

    }

    for(int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        //left boundary condition
        discretization_->v(discretization_->vIBegin(),j) = 2.0*settings_.dirichletBcLeft[1]-discretization_->v(discretization_->vIBegin()+1,j);

        //right boundary condition
        discretization_->v(discretization_->vIEnd()-1,j)  = 2.0*settings_.dirichletBcRight[1]-discretization_->v(discretization_->vIEnd()-2,j);
    }
   
}

void Computation::applyBoundaryValuesFandG(){

    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        discretization_->f(i, discretization_->uJBegin()) = discretization_->u(i, discretization_->uJBegin());
        discretization_->f(i, discretization_->uJEnd() - 1) = discretization_->u(i, discretization_->uJEnd() - 1);
    }


    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        discretization_->g(i, discretization_->vJBegin()) = discretization_->v(i, discretization_->vJBegin());
        discretization_->g(i, discretization_->vJEnd() - 1) = discretization_->v(i, discretization_->vJEnd() - 1);
    }

    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
        discretization_->f(discretization_->uIBegin(), j) = discretization_->u(discretization_->uIBegin(), j);
        discretization_->f(discretization_->uIEnd() - 1, j) = discretization_->u(discretization_->uIEnd() - 1, j);
    }

    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
        discretization_->g(discretization_->vIBegin(), j) = discretization_->v(discretization_->vIBegin(), j);
        discretization_->g(discretization_->vIEnd() - 1, j) = discretization_->v(discretization_->vIEnd() - 1, j);
    }
}

void Computation::computePreliminaryVelocities()
{

    for (int i = discretization_->uIBegin()+1; i < discretization_->uIEnd()-1; i++)
    {
        for (int j = discretization_->uJBegin()+1; j < discretization_->uJEnd()-1; j++)
        {   
            double diffusion_u       = discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j);
            double convection_u      = discretization_->computeDu2Dx(i,j)  + discretization_->computeDuvDy(i,j);
            
            discretization_->f(i,j) = discretization_->u(i,j) + dt_*(diffusion_u/settings_.re - convection_u + settings_.g[0]);
        }
        
    }

    for (int i = discretization_->vIBegin()+1; i < discretization_-> vIEnd()-1; i++)
    {
        for (int j = discretization_->vJBegin()+1; j < discretization_-> vJEnd()-1; j++)
        {
            double diffusion_v  = discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j);
            double convection_v = discretization_->computeDv2Dy(i,j)  + discretization_->computeDuvDx(i,j);
            
            discretization_->g(i,j) = discretization_->v(i,j) + dt_*(diffusion_v/settings_.re - convection_v + settings_.g[1]);

        }
    }
}

void Computation::computeRightHandSide()
{
    for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++)
    {
        for (int j = discretization_->pJBegin()+1; j < discretization_-> pJEnd()-1; j++)
        {
            double dfdx = (discretization_->f(i,j) - discretization_-> f(i-1,j))/discretization_->dx();
            double dgdy = (discretization_->g(i,j) - discretization_-> g(i,j-1))/discretization_->dy();

            discretization_->rhs(i,j) = (dfdx + dgdy)/dt_;
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
    for(int i = discretization_->uIBegin()+1; i < discretization_->uIEnd()-1; i++){
        for(int j = discretization_->uJBegin()+1; j < discretization_->uJEnd()-1; j++){
            discretization_->u(i,j) = discretization_->f(i,j)-dt_*discretization_->computeDpDx(i,j);
        }
    }

    //set values for v
    for(int i = discretization_->vIBegin()+1; i < discretization_->vIEnd()-1; i++){
        for(int j = discretization_->vJBegin()+1; j < discretization_-> vJEnd()-1; j++){
            discretization_->v(i,j) = discretization_->g(i,j)-dt_*discretization_->computeDpDy(i,j);
        }
    }
}

