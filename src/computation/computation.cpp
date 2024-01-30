#include "computation/computation.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/sor.h"
#include "discretization/central_grid.h"

/**
 * Initialize the computation object
 * 
 * Parse the settings from the parameter file that is given as the command line argument
 * It implements the time stepping scheme, computes all the terms and calls the pressure solver.
 */
// Computation::Computation(): ci_(settings_.Q),wi_(settings_.Q)
// {

// }
void Computation::initialize(std::string filename)
{
    settings_ = Settings();
    // Load settings from file
  std::cout <<"#########################"<<std::endl;    
    //settings_.loadFromFile(filename);
    // Print settings

    //settings_.printSettings();
 

    // Initialize discretization
    for (int i = 0; i < 2; i++)
        meshWidth_[i] = settings_.physicalSize[i] / settings_.nCells[i];

    // if (settings_.useDonorCell) {
    //     discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    // }
    // else {
    //     discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    // }

    // // Initialize solver
    // if (settings_.pressureSolver == "SOR") {
    //     pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon,
    //                                             settings_.maximumNumberOfIterations, settings_.omega);
    // }
    // else if (settings_.pressureSolver == "GaussSeidel") {
    //     pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon,
    //                                             settings_.maximumNumberOfIterations);

    // } else {
    //     std::cout << "Solver not found!" << std::endl;
    // }

    //if (settings_.useLBM){
        doLBM = true;
        cdiscretization_ = std::make_shared<LbmDiscretization>(settings_.nCells, meshWidth_);
        cs_ =(double) 1./sqrt(3);
        nu_ = settings_.dirichletBcTop[0]*settings_.L_lbm/settings_.re;
        tau_ = 1.;//0.5+nu_/(cs_*cs_);
        dt_ =1.;//settings_.re*nu_/(settings_.L_lbm*settings_.L_lbm);
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);        // Initialize output writers
        outputWriterText_ = std::make_unique<OutputWriterText>(cdiscretization_);
        outputWriterParaview_ = std::make_unique<OutputWriterParaview>(cdiscretization_);    
        Collision();
        for(int i=cdiscretization_->pIBegin(); i < cdiscretization_->pIEnd();i++)
            {
                for(int j=cdiscretization_->pJBegin(); j < cdiscretization_->pJEnd();j++)
                    {
                        cdiscretization_->rho(i,j) = 1.;
                        for(int k = 0; k<9; k++){
                            cdiscretization_->pdf(i,j,k)  = cdiscretization_->pdfeq(i,j,k);      
                        }
                    }
            }     

    // }
    // else 
    // {
    //     doLBM = false;
    //     // Initialize output writers
    //     outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
    //     outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);        
    // }

};

/**
 * Run the whole simulation until tend
 */
void Computation::runSimulation() {
    int t_iter = 0;
    double time = 0.0;
    while (time < settings_.endTime){

//         if(doLBM == false){
//         /*
//         * 1) Apply boundary values (for u, v, F, G)
//         */
//         applyBoundaryValues();
//         applyPreliminaryBoundaryValues();

//         /*
//         * 2) Compute the next time step width
//         */
//         computeTimeStepWidth();
//         // endTime should be reached exactly:
//         if (time + dt_ > settings_.endTime) {
//             dt_ = settings_.endTime - time;
//         }
//         time += dt_;

//         /*
//         * 3) Compute preliminary velocities (F, G)
//         */
//         computePreliminaryVelocities();

//         /*
//         * 4) Compute the right hand side (rhs)
//         */
//         computeRightHandSide();

//         /*
//         * 5) Compute the pressure (p) by solving the Poisson equation
//         */
//         computePressure();

//         /*
//         * 6) Update the velocities (u, v)
//         */
//         computeVelocities();

//         /*
//         * 7) Output debug information and simulation results
//         */
// #ifndef NDEBUG
//         cout << "time step " << t_iter << ", t: " << time << "/" << settings_.endTime << ", dt: " << dt_ <<
//             ", res. " << pressureSolver_->residualNorm() << ", solver iterations: " << pressureSolver_->iterations() << endl;
//         //outputWriterText_->writePressureFile();
//         outputWriterText_->writeFile(time);
// #endif
//         outputWriterParaview_->writeFile(time);
//         }
//         else
//         { 



        time += dt_;




        settings_.dirichletBcTop[0] = 1000.*(1.0-exp(-t_iter*t_iter/(2.0*10*settings_.nCells[0]*10*settings_.nCells[0])));
                t_iter++;
        Collision();
        FillVelocitiesAndPressure();
        //std::cout<<"+++++++++++++++++++"<<std::endl;     

        Streaming();
        FillVelocitiesAndPressure();
        LBMapplyBoundaryValues();        
        //FillVelocitiesAndPressure();         
        //std::cout<<"*******************"<<std::endl;

        //LBMBounceBack();

        //std::cout<<"-------------------"<<std::endl;      
        FillVelocitiesAndPressure();        
        //std::cout<<"+++++++++++++++++++"<<std::endl;
        outputWriterText_->writeFile(time);        
        outputWriterParaview_->writeFile(time);


        // std::cout<<"                 "<<std::endl;
        // for(int k = 0; k<9;k++){
        //     std::cout<<"  dt:  "<<dt_<<"  "<<"########pdf0/2: "<<cdiscretization_->pdf(0,4,k)<<" "<<"           ########pdf0/1: "<<cdiscretization_->pdf(1,2,k)<<" "<<"            ########pdf0/0: "<<cdiscretization_->pdf(1,3,k)<<std::endl;
        // }



        //}
    }

};

void Computation::LBMapplyBoundaryValues() {
    //top condition
    int pdfJEnd = cdiscretization_->pdfJEnd()-1;
    int pdfJBegin = cdiscretization_->pdfJBegin();
    int pdfIBegin = cdiscretization_->pdfIBegin();
    int pdfIEnd = cdiscretization_->pdfIEnd()-1;
//top
    for(int i = cdiscretization_->pdfIBegin(); i < cdiscretization_->pdfIEnd(); i++) {

        double rho_N = 1./(1.+settings_.dirichletBcTop[1])*(cdiscretization_->pdf(i,pdfJEnd,0)+cdiscretization_->pdf(i,pdfJEnd,1)+
            cdiscretization_->pdf(i,pdfJEnd,3)+2.*(cdiscretization_->pdf(i,pdfJEnd,2)+cdiscretization_->pdf(i,pdfJEnd,6)+cdiscretization_->pdf(i,pdfJEnd,5)));
        cdiscretization_->rho(i,pdfJEnd) = rho_N; 
        cdiscretization_->pdf(i,pdfJEnd,4) = cdiscretization_->pdf(i,pdfJEnd,2)-2./3.*rho_N*cdiscretization_->u(i,pdfJEnd);
        cdiscretization_->pdf(i,pdfJEnd,7) = (cdiscretization_->pdf(i,pdfJEnd,5)+1./2.*(cdiscretization_->pdf(i,pdfJEnd,1)-cdiscretization_->pdf(i,pdfJEnd,3)) 
                                                -1./6.*rho_N*cdiscretization_->v(i,pdfJEnd)-1./2.*rho_N*settings_.dirichletBcTop[0]);
        cdiscretization_->pdf(i,pdfJEnd,8) = (cdiscretization_->pdf(i,pdfJEnd,6)-1./2.*(cdiscretization_->pdf(i,pdfJEnd,1)-cdiscretization_->pdf(i,pdfJEnd,3))
                                                -1./6.*rho_N*cdiscretization_->v(i,pdfJEnd)+1./2.*rho_N*settings_.dirichletBcTop[0]);

    }
//Bottom
    for(int i = cdiscretization_->pdfIBegin(); i < cdiscretization_->pdfIEnd(); i++) {

        double rho_N = 1./(1.-settings_.dirichletBcBottom[1])*(cdiscretization_->pdf(i,pdfJBegin,0)+cdiscretization_->pdf(i,pdfJBegin,1)+
            cdiscretization_->pdf(i,pdfJBegin,3)+2.*(cdiscretization_->pdf(i,pdfJBegin,4)+cdiscretization_->pdf(i,pdfJBegin,7)+cdiscretization_->pdf(i,pdfJBegin,8)));
        cdiscretization_->rho(i,pdfJBegin) = rho_N; 
        cdiscretization_->pdf(i,pdfJBegin,2) = cdiscretization_->pdf(i,pdfJBegin,4)+2./3.*rho_N*settings_.dirichletBcBottom[1];
        cdiscretization_->pdf(i,pdfJBegin,5) = (cdiscretization_->pdf(i,pdfJBegin,7)-1./2.*(cdiscretization_->pdf(i,pdfJBegin,1)-cdiscretization_->pdf(i,pdfJBegin,3)) 
                                                +1./6.*rho_N*settings_.dirichletBcBottom[1]+1./2.*rho_N*settings_.dirichletBcBottom[0]);
        cdiscretization_->pdf(i,pdfJBegin,6) = (cdiscretization_->pdf(i,pdfJBegin,8)+1./2.*(cdiscretization_->pdf(i,pdfJBegin,1)-cdiscretization_->pdf(i,pdfJBegin,3))
                                                +1./6.*rho_N*settings_.dirichletBcBottom[1]-1./2.*rho_N*settings_.dirichletBcBottom[0]);

    }
//left
    for(int j = cdiscretization_->pdfJBegin(); j < cdiscretization_->pdfJEnd(); j++) {

        double rho_N = 1./(1.-settings_.dirichletBcLeft[0])*(cdiscretization_->pdf(pdfIBegin,j,0)+cdiscretization_->pdf(pdfIBegin,j,2)+
            cdiscretization_->pdf(pdfIBegin,j,4)+2.*(cdiscretization_->pdf(pdfIBegin,j,3)+cdiscretization_->pdf(pdfIBegin,j,6)+cdiscretization_->pdf(pdfIBegin,j,7)));
        cdiscretization_->rho(pdfIBegin,j) = rho_N; 
        cdiscretization_->pdf(pdfIBegin,j,1) = cdiscretization_->pdf(pdfIBegin,j,3)+2./3.*rho_N*settings_.dirichletBcLeft[0];
        cdiscretization_->pdf(pdfIBegin,j,5) = (cdiscretization_->pdf(pdfIBegin,j,7)-1./2.*(cdiscretization_->pdf(pdfIBegin,j,2)-cdiscretization_->pdf(pdfIBegin,j,4)) 
                                                +1./6.*rho_N*settings_.dirichletBcLeft[0]+1./2.*rho_N*settings_.dirichletBcLeft[1]);
        cdiscretization_->pdf(pdfIBegin,j,8) = (cdiscretization_->pdf(pdfIBegin,j,6)+1./2.*(cdiscretization_->pdf(pdfIBegin,j,2)-cdiscretization_->pdf(pdfIBegin,j,4))
                                                +1./6.*rho_N*settings_.dirichletBcLeft[0]-1./2.*rho_N*settings_.dirichletBcLeft[1]);

    }
//right
    for(int j = cdiscretization_->pdfJBegin(); j < cdiscretization_->pdfJEnd(); j++) {

        double rho_N = 1./(1.+settings_.dirichletBcRight[0])*(cdiscretization_->pdf(pdfIEnd,j,0)+cdiscretization_->pdf(pdfIEnd,j,2)+
            cdiscretization_->pdf(pdfIEnd,j,4)+2.*(cdiscretization_->pdf(pdfIEnd,j,5)+cdiscretization_->pdf(pdfIEnd,j,1)+cdiscretization_->pdf(pdfIEnd,j,8)));
        cdiscretization_->rho(pdfIEnd,j) = rho_N; 
        cdiscretization_->pdf(pdfIEnd,j,3) = cdiscretization_->pdf(pdfIEnd,j,1)-2./3.*rho_N*settings_.dirichletBcRight[0];
        cdiscretization_->pdf(pdfIEnd,j,7) = (cdiscretization_->pdf(pdfIEnd,j,5)+1./2.*(cdiscretization_->pdf(pdfIEnd,j,2)-cdiscretization_->pdf(pdfIEnd,j,4)) 
                                                -1./6.*rho_N*settings_.dirichletBcRight[0]-1./2.*rho_N*settings_.dirichletBcRight[1]);
        cdiscretization_->pdf(pdfIEnd,j,6) = (cdiscretization_->pdf(pdfIEnd,j,8)-1./2.*(cdiscretization_->pdf(pdfIEnd,j,2)-cdiscretization_->pdf(pdfIEnd,j,4))
                                                -1./6.*rho_N*settings_.dirichletBcRight[0]+1./2.*rho_N*settings_.dirichletBcRight[1]);

    }
//set corners
//bottom left
cdiscretization_->u(0,0) = cdiscretization_->u(1,0);
cdiscretization_->v(0,0) = cdiscretization_->v(1,0);
cdiscretization_->rho(0,0) = cdiscretization_->rho(1,0);
cdiscretization_->pdf(0,0,1) = cdiscretization_->pdf(0,0,3)+(2.0/3.0)*cdiscretization_->rho(0,0)*cdiscretization_->u(1,0);
cdiscretization_->pdf(0,0,2) = cdiscretization_->pdf(0,0,4)+(2.0/3.0)*cdiscretization_->rho(0,0)*cdiscretization_->v(1,0);
cdiscretization_->pdf(0,0,5) = cdiscretization_->pdf(0,0,7)+(1.0/6.0)*cdiscretization_->rho(0,0)*cdiscretization_->u(1,0)+(1./6.)*cdiscretization_->rho(0,0)*cdiscretization_->v(1,0);
cdiscretization_->pdf(0,0,6) = 0.0;
cdiscretization_->pdf(0,0,8) = 0.0;
cdiscretization_->pdf(0,0,0) = cdiscretization_->rho(0,0)-cdiscretization_->pdf(0,0,1)-cdiscretization_->pdf(0,0,2)-cdiscretization_->pdf(0,0,3)-cdiscretization_->pdf(0,0,4)
                                -cdiscretization_->pdf(0,0,5)-cdiscretization_->pdf(0,0,6)-cdiscretization_->pdf(0,0,7)-cdiscretization_->pdf(0,0,8);
//top left
cdiscretization_->u(0,pdfJEnd-1) = cdiscretization_->u(1,pdfJEnd-1);
cdiscretization_->v(0,pdfJEnd-1) = cdiscretization_->v(1,pdfJEnd-1);
cdiscretization_->rho(0,cdiscretization_->pdfJEnd()-2) = cdiscretization_->rho(1,cdiscretization_->pdfJEnd()-2);
cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-2,1) = cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-2,3)
                                                            +(2.0/3.0)*cdiscretization_->rho(0,cdiscretization_->pdfJEnd()-1)*cdiscretization_->u(1,pdfJEnd-1);
cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,4) = cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,2)
                                                            -(2.0/3.0)*cdiscretization_->rho(0,cdiscretization_->pdfJEnd()-1)*cdiscretization_->v(1,pdfJEnd-1);
cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,8) = cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,6)
                                                            +(1.0/6.0)*cdiscretization_->rho(0,cdiscretization_->pdfJEnd()-1)*cdiscretization_->u(1,pdfJEnd-1)
                                                            -(1./6.)*cdiscretization_->rho(0,cdiscretization_->pdfJEnd()-1)*cdiscretization_->v(1,pdfJEnd-1);
cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,5) = 0.0;
cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,7) = 0.0;
cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,0) = cdiscretization_->rho(0,cdiscretization_->pdfJEnd()-1)-cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,1)
                                -cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,2)-cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,3)-cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,4)
                                -cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,5)-cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,6)
                                -cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,7)-cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,8);
//top right
cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1) = cdiscretization_->rho(cdiscretization_->pdfIEnd()-2,cdiscretization_->pdfJEnd()-2);
cdiscretization_->u(pdfIEnd,pdfJEnd) = cdiscretization_->u(pdfIEnd-1,pdfJEnd-1);
cdiscretization_->v(pdfIEnd,pdfJEnd) = cdiscretization_->v(pdfIEnd-1,pdfJEnd-1);

cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,3) = cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,1)
                                                            -(2.0/3.0)*cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1)*cdiscretization_->u(pdfIEnd,pdfJEnd);
cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,4) = cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,2)
                                                            -(2.0/3.0)*cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1)*cdiscretization_->v(pdfIEnd,pdfJEnd);
cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,7) = cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,5)
                                                            -(1.0/6.0)*cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1)*cdiscretization_->u(pdfIEnd,pdfJEnd)
                                                            -(1./6.)*cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1)*cdiscretization_->v(pdfIEnd,pdfJEnd);
cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,6) = 0.0;
cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,8) = 0.0;
cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,0) = cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1)
                                                                                        -cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,1)
                                -cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,2)
                                -cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,3)-cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,4)
                                -cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,5)-cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,6)
                                -cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,7)-cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,8);
//bottom right
cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,0) = cdiscretization_->rho(cdiscretization_->pdfIEnd()-2,0);
cdiscretization_->u(pdfIEnd,pdfJBegin) = cdiscretization_->u(pdfIEnd-1,pdfJBegin+1);
cdiscretization_->v(pdfIEnd,pdfJBegin) = cdiscretization_->v(pdfIEnd-1,pdfJBegin+1);
cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,3) = cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,1)
                                                            -(2.0/3.0)*cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,0)*cdiscretization_->u(pdfIEnd,pdfJBegin);
cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,2) = cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,4)
                                                            +(2.0/3.0)*cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,0)*cdiscretization_->v(pdfIEnd,pdfJBegin);
cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,6) = cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,8)
                                                            -(1.0/6.0)*cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,0)*cdiscretization_->u(pdfIEnd,pdfJBegin)
                                                            +(1./6.)*cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,0)*cdiscretization_->v(pdfIEnd,pdfJBegin);
cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,5) = 0.0;
cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,7) = 0.0;
cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,0) = cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,0)-cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,1)
                                                            -cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,2)-cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,3)-cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,4)
                                                            -cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,5)-cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,6)
                                                            -cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,7)-cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,8);
};

void Computation::LBMBounceBack() {
//    top condition
    int pdfJEnd = cdiscretization_->pdfJEnd()-1;
    int pdfJBegin = cdiscretization_->pdfJBegin();
    int pdfIBegin = cdiscretization_->pdfIBegin();
    int pdfIEnd = cdiscretization_->pdfIEnd()-1;

    for(int i = cdiscretization_->pdfIBegin(); i < cdiscretization_->pdfIEnd(); i++) {
        // cdiscretization_->pdf(i,pdfJEnd,4) =  cdiscretization_->pdf(i,pdfJEnd,2);
        // cdiscretization_->pdf(i,pdfJEnd,7) =  (cdiscretization_->pdf(i,pdfJEnd,5));
        // cdiscretization_->pdf(i,pdfJEnd,8) =  (cdiscretization_->pdf(i,pdfJEnd,6));

        cdiscretization_->pdf(i,pdfJBegin,2) =  cdiscretization_->pdf(i,pdfJBegin,4);
        cdiscretization_->pdf(i,pdfJBegin,5) =  (cdiscretization_->pdf(i,pdfJBegin,7));
        cdiscretization_->pdf(i,pdfJBegin,6) =  (cdiscretization_->pdf(i,pdfJBegin,8));
    }
    for(int j = cdiscretization_->pdfIBegin(); j < cdiscretization_->pdfIEnd(); j++)
    {
        cdiscretization_->pdf(pdfIBegin,j,1) = cdiscretization_->pdf(pdfIBegin,j,3); 
        cdiscretization_->pdf(pdfIBegin,j,8) = cdiscretization_->pdf(pdfIBegin,j,6);
        cdiscretization_->pdf(pdfIBegin,j,5) = cdiscretization_->pdf(pdfIBegin,j,7);

        cdiscretization_->pdf(pdfIEnd,j,3) = cdiscretization_->pdf(pdfIEnd,j,1);
        cdiscretization_->pdf(pdfIEnd,j,6) = cdiscretization_->pdf(pdfIEnd,j,8);
        cdiscretization_->pdf(pdfIEnd,j,7) = cdiscretization_->pdf(pdfIEnd,j,5);         
    }

};


void Computation::FillVelocitiesAndPressure(){
    for(int i=cdiscretization_->pIBegin(); i < cdiscretization_->pIEnd();i++)
    {
        for(int j=cdiscretization_->pJBegin(); j < cdiscretization_->pJEnd();j++)
        {
              double p = 0.;
              double u = 0.;
              double v = 0.;
              double rho = 0.;
            for(int k =0;k<9;k++)
            {
              rho += cdiscretization_->pdf(i,j,k);
              p += cdiscretization_->pdf(i,j,k)*cs_*cs_;                
              u += cdiscretization_->pdf(i,j,k)*cdiscretization_->ci_x(k);
              v += cdiscretization_->pdf(i,j,k)*cdiscretization_->ci_y(k);  
            //std::cout<<"###### "<<cdiscretization_->pdf(i,j,k)<< "####rho "<<cdiscretization_->rho(i,j)<<std::endl;               
            }

             cdiscretization_->rho(i,j) = rho;
             cdiscretization_->p(i,j) = p;
             cdiscretization_->u(i,j) = u/rho;
             cdiscretization_->v(i,j) = v/rho;
        }
    }
};

void Computation::Collision(){
    int pdfJEnd = cdiscretization_->pdfJEnd();
    int pdfJBegin = cdiscretization_->pdfJBegin();
    int pdfIBegin = cdiscretization_->pdfIBegin();
    int pdfIEnd = cdiscretization_->pdfIEnd();

        for(int i = pdfIBegin; i < pdfIEnd; i++){
            for(int j = pdfJBegin; j<pdfJEnd; j++){
                for(int k=0;k<9;k++){
                    cdiscretization_->pdfeq(i,j,k) = cdiscretization_->wi(k)*cdiscretization_->rho(i,j)*(1.+(cdiscretization_->u(i,j)*cdiscretization_->ci_x(k)
                                                                                    +cdiscretization_->v(i,j)*cdiscretization_->ci_y(k))/(cs_*cs_)
                                                        +1./2.*((cdiscretization_->u(i,j)*cdiscretization_->ci_x(k)
                                                                +cdiscretization_->v(i,j)*cdiscretization_->ci_y(k))/(cs_*cs_))
                                                        *((cdiscretization_->u(i,j)*cdiscretization_->ci_x(k)
                                                                +cdiscretization_->v(i,j)*cdiscretization_->ci_y(k))/(cs_*cs_))
                                                        -1./2.*(cdiscretization_->u(i,j)*cdiscretization_->u(i,j)
                                                                +cdiscretization_->v(i,j)*cdiscretization_->v(i,j))/(cs_*cs_));
                //std::cout<<"###### "<<cdiscretization_->pdfeq(i,j,k)<< "####rho "<<cdiscretization_->rho(i,j)<<std::endl;
                cdiscretization_->pdf(i,j,k) = cdiscretization_->pdf(i,j,k)-dt_/tau_*(cdiscretization_->pdf(i,j,k)-cdiscretization_->pdfeq(i,j,k));
                }
            }
        }

};

void Computation::Streaming(){
    std::array<int, 9> cx_ = {0,1,0,-1,0,1,-1,-1,1};
    std::array<int, 9> cy_ = {0,0,1,0,-1,1,1,-1,-1}; 
    //write pdf in pdfold
    for(int i=cdiscretization_->pdfIBegin(); i < cdiscretization_->pdfIEnd();i++)
    {
        for(int j=cdiscretization_->pdfJBegin(); j < cdiscretization_->pdfJEnd();j++)
        {
            for(int k = 0;k<9;k++){
            cdiscretization_->pdfold(i,j,k) = cdiscretization_->pdf(i,j,k);
            }
        }
    }
    int pdfIBeginInner = cdiscretization_->pdfIBegin()+1;
    int pdfIEndInner = cdiscretization_->pdfIEnd()-1;
    int pdfJBeginInner = cdiscretization_->pdfJBegin()+1;
    int pdfJEndInner = cdiscretization_->pdfJEnd()-1;

   for(int i = pdfIBeginInner-1; i<pdfIEndInner+1;i++){
        for(int j = pdfJBeginInner-1; j<pdfJEndInner+1;j++){
                    cdiscretization_->pdf(i,j,0) = cdiscretization_->pdfold(i,j,0);
        }
    }
    for(int i = pdfIBeginInner-1; i<pdfIEndInner;i++){
        for(int j = pdfJBeginInner-1; j<pdfJEndInner+1;j++){
                    cdiscretization_->pdf(i+1,j,1) = cdiscretization_->pdfold(i,j,1);
        }
    }
    for(int i = pdfIBeginInner; i<pdfIEndInner+1;i++){
        for(int j = pdfJBeginInner-1; j<pdfJEndInner+1;j++){
                    cdiscretization_->pdf(i-1,j,3) = cdiscretization_->pdfold(i,j,3);
        }
    }
    for(int i = pdfIBeginInner-1; i<pdfIEndInner+1;i++){
        for(int j = pdfJBeginInner-1; j<pdfJEndInner;j++){
                    cdiscretization_->pdf(i,j+1,2) = cdiscretization_->pdfold(i,j,2);
        }
    }
    for(int i = pdfIBeginInner-1; i<pdfIEndInner+1;i++){
        for(int j = pdfJBeginInner; j<pdfJEndInner+1;j++){
                    cdiscretization_->pdf(i,j-1,4) = cdiscretization_->pdfold(i,j,4);
        }
    }
    for(int i = pdfIBeginInner-1; i<pdfIEndInner;i++){
        for(int j = pdfJBeginInner-1; j<pdfJEndInner;j++){
                    cdiscretization_->pdf(i+1,j+1,5) = cdiscretization_->pdfold(i,j,5);
        }
    }
    for(int i = pdfIBeginInner; i<pdfIEndInner+1;i++){
        for(int j = pdfJBeginInner; j<pdfJEndInner+1;j++){
                    cdiscretization_->pdf(i-1,j-1,7) = cdiscretization_->pdfold(i,j,7);
        }
    }
    for(int i = pdfIBeginInner; i<pdfIEndInner+1;i++){
        for(int j = pdfJBeginInner-1; j<pdfJEndInner;j++){
                    cdiscretization_->pdf(i-1,j+1,6) = cdiscretization_->pdfold(i,j,6);
        }
    }
    for(int i = pdfIBeginInner-1; i<pdfIEndInner;i++){
        for(int j = pdfJBeginInner; j<pdfJEndInner+1;j++){
                    cdiscretization_->pdf(i+1,j-1,8) = cdiscretization_->pdfold(i,j,8);
        }
    }     
};

/**
 * Set the boundary values of the velocities (u, v)
 * 
 * Left and right boundaries should overwrite bottom and top boundaries
 */
void Computation::applyBoundaryValues() {
    // set boundary values for u at bottom and top side (lower priority)
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        // set boundary values for u at bottom side
        discretization_->u(i, discretization_->uJBegin()) =
                2.0 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uInteriorJBegin());
        // set boundary values for u at top side
        discretization_->u(i, discretization_->uJEnd() - 1) =
                2.0 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uInteriorJEnd() - 1);
    }

    // set boundary values for v at bottom and top side (lower priority)
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        // set boundary values for v at bottom side
        discretization_->v(i, discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
        // set boundary values for v at top side
        discretization_->v(i, discretization_->vJEnd() - 1) = settings_.dirichletBcTop[1];
    }

    // set boundary values for u at left and right side (higher priority)
    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
        // set boundary values for u at left side
        discretization_->u(discretization_->uIBegin(), j) = settings_.dirichletBcLeft[0];
        // set boundary values for u at right side
        discretization_->u(discretization_->uIEnd() - 1, j) = settings_.dirichletBcRight[0];
    }

    // set boundary values for v at left and right side (higher priority)
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
        // set boundary values for v at left side
        discretization_->v(discretization_->vIBegin(), j) =
                2.0 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin() + 1, j);
        // set boundary values for v at right side
        discretization_->v(discretization_->vIEnd() - 1, j) =
                2.0 * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vInteriorIEnd() - 1, j);
    }
};

/**
 * Set the boundary values of the preliminary velocities (u, v)
 * 
 * Left and right boundaries should overwrite bottom and top boundaries
 */
void Computation::applyPreliminaryBoundaryValues(){
    // set boundary values for F at bottom and top side (lower priority)
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
        // set boundary values for F at bottom side
        discretization_->f(i, discretization_->uJBegin()) = discretization_->u(i, discretization_->uJBegin());
        // set boundary values for F at top side
        discretization_->f(i, discretization_->uJEnd() - 1) = discretization_->u(i, discretization_->uJEnd() - 1);
    }

    // set boundary values for G at bottom and top side (lower priority)
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
        // set boundary values for v at bottom side
        discretization_->g(i, discretization_->vJBegin()) = discretization_->v(i, discretization_->vJBegin());
        // set boundary values for v at top side
        discretization_->g(i, discretization_->vJEnd() - 1) = discretization_->v(i, discretization_->vJEnd() - 1);
    }

    // set boundary values for F at left and right side (higher priority)
    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++) {
        // set boundary values for F at left side
        discretization_->f(discretization_->uIBegin(), j) = discretization_->u(discretization_->uIBegin(), j);
        // set boundary values for F at right side
        discretization_->f(discretization_->uIEnd() - 1, j) = discretization_->u(discretization_->uIEnd() - 1, j);
    }

    // set boundary values for G at left and right side (higher priority)
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++) {
        // set boundary values for G at left side
        discretization_->g(discretization_->vIBegin(), j) = discretization_->v(discretization_->vIBegin(), j);
        // set boundary values for G at right side
        discretization_->g(discretization_->vIEnd() - 1, j) = discretization_->v(discretization_->vIEnd() - 1, j);
    }
};

/**
 * Compute the preliminary velocities (F, G) using finite differences
 */ 
void Computation::computePreliminaryVelocities(){
    // Compute F in the interior of the domain
    for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            double lap_u = discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j);
            double conv_u = discretization_->computeDu2Dx(i,j) + discretization_->computeDuvDy(i,j);
            discretization_->f(i,j) = discretization_->u(i,j) + dt_ * (lap_u / settings_.re - conv_u + settings_.g[0]);
        }
    }

    // Compute G in the interior of the domain
    for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            double lap_v = discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j);
            double conv_v = discretization_->computeDv2Dy(i,j) + discretization_->computeDuvDx(i,j);
            discretization_->g(i,j) = discretization_->v(i,j) + dt_ * (lap_v / settings_.re - conv_v + settings_.g[1]);
        }
    }
};

/**
 * Compute the pressure p by solving the Poisson equation
 */
void Computation::computePressure() {
    pressureSolver_->solve();
};

/**
 * Compute the right hand side rhs of the pressure Poisson equation 
 */
void Computation::computeRightHandSide() {
    // Compute rhs in the interior of the domain using finite differences
    for (int i = discretization_->rhsInteriorIBegin(); i < discretization_->rhsInteriorIEnd(); i++) {
        for (int j = discretization_->rhsInteriorJBegin(); j < discretization_->rhsInteriorJEnd(); j++) {
            double fx = (discretization_->f(i, j) - discretization_->f(i - 1, j)) / discretization_->dx();
            double gy = (discretization_->g(i, j) - discretization_->g(i, j - 1)) / discretization_->dy();
            discretization_->rhs(i, j) = (fx + gy) / dt_;
        }
    }
};

/**
 * Compute the time step width dt based on the maximum velocities
 */
void Computation::computeTimeStepWidth() {
    // Compute maximal time step width regarding the diffusion
    double dt_diff = settings_.re / 2 / (1 / (discretization_->dx() * discretization_->dx()) + 1 / (discretization_->dy() * discretization_->dy()) );

    // Compute maximal time step width regarding the convection u
    double dt_conv_u = discretization_->dx() / discretization_->u().absMax();
    
    // Compute maximal time step width regarding the convection v
    double dt_conv_v = discretization_->dy() / discretization_->v().absMax();

    // Set the appropriate time step width by using a security factor tau
    dt_ = settings_.tau * std::min({dt_diff, dt_conv_u, dt_conv_v});
};

/**
 * Compute the new velocities (u, v) based on the preliminary velocities (F, G) and the pressure (p)
 */
void Computation::computeVelocities() {
    // Compute u in the interior of the domain
    for (int i = discretization_->uInteriorIBegin(); i < discretization_->uInteriorIEnd(); i++) {
        for (int j = discretization_->uInteriorJBegin(); j < discretization_->uInteriorJEnd(); j++) {
            discretization_->u(i, j) = discretization_->f(i, j) - dt_ * discretization_->computeDpDx(i, j);
        }
    }

    // Compute v in the interior of the domain
    for (int i = discretization_->vInteriorIBegin(); i < discretization_->vInteriorIEnd(); i++) {
        for (int j = discretization_->vInteriorJBegin(); j < discretization_->vInteriorJEnd(); j++) {
            discretization_->v(i, j) = discretization_->g(i, j) - dt_ * discretization_->computeDpDy(i, j);
        }
    }
};
