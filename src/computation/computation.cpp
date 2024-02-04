#include "computation/computation.h"
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
    // settings_.loadFromFile(filename);
    // Print settings

    // settings_.printSettings();

    // Initialize discretization
    for (int i = 0; i < 3; i++)
        meshWidth_[i] = settings_.physicalSize[i] / settings_.nCells[i];

    // doLBM = true;
    cdiscretization_ = std::make_shared<LbmDiscretization>(settings_.nCells, meshWidth_);

    cs_ = (double)1. / sqrt(3);
    nu_ = settings_.dirichletBcTop[0] * settings_.L_lbm / settings_.re;
    tau_ = 0.5 + nu_ / (cs_ * cs_);
    // tau_ =0.65;
    dt_ = 1.; // settings_.re*nu_/(settings_.nCells[0]*settings_.nCells[0]);
    std::cout << dt_ << std::endl;
    //(settings_.nCells[1]/settings_.re)*settings_.dirichletBcTop[0]+.5;

    // tau_m = 1./4.*1/(tau_-0.5)+0.5;
    //  tau_om_m = 1./tau_m;
    //  tau_om_p = 1./tau_;
    //  tau_om   = 1./tau_;
    //  dt_ = meshWidth_[0]*settings_.dirichletBcTop[0];
    //  tau_ = 1.25;

    std::cout << tau_ << std::endl;
    outputWriterText_ = std::make_unique<OutputWriterText>(cdiscretization_);
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(cdiscretization_);
    // double v = meshWidth_[0]/dt_;

    // for(int j=0; j < 9;j++)
    //     {
    //         double holderx = 0.0;
    //         double holdery = 0.0;
    //         holderx = cdiscretization_->ci_x(j)*v;
    //         holdery = cdiscretization_->ci_y(j)*v;
    //         cdiscretization_->ci_x(j) = holderx;
    //         cdiscretization_->ci_y(j) =holdery;
    //     }
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
void Computation::runSimulation()
{
    int t_iter = 0;
    double time = 0.0;
    while (time < settings_.endTime)
    {

        // settings_.dirichletBcTop[0] = 1000.*(1.0-exp(-t_iter*t_iter/(2.0*10*settings_.nCells[0]*10*settings_.nCells[0])));
        //         t_iter++;

        // std::cout<<"+++++++++++++++++++"<<std::endl;
        if (time == 0.0)
        {
            for (int i = cdiscretization_->pIBegin(); i < cdiscretization_->pIEnd(); i++)
            {
                for(int j = cdiscretization_->pIBegin(); j <cdiscretization_->pIEnd(); j++){
                cdiscretization_->u(i,j, cdiscretization_->uJBegin()) = settings_.dirichletBcBottom[0];
                cdiscretization_->u(i,j,cdiscretization_->uJEnd() - 1) = settings_.dirichletBcTop[0];
                cdiscretization_->v(i,j, cdiscretization_->vJBegin()) = settings_.dirichletBcBottom[1];
                cdiscretization_->v(i,j, cdiscretization_->vJEnd() - 1) = settings_.dirichletBcTop[1];
                cdiscretization_->v(i,j, cdiscretization_->vJBegin()) = settings_.dirichletBcBottom[1];
                cdiscretization_->v(i,j, cdiscretization_->vJEnd() - 1) = settings_.dirichletBcTop[1];
                }
            }
            for (int j = cdiscretization_->pJBegin(); j < cdiscretization_->pJEnd(); j++)
            {
                //cdiscretization_->u(cdiscretization_->uJBegin(), j) = settings_.dirichletBcLeft[0];
                //cdiscretization_->u(cdiscretization_->uJEnd() - 1, j) = settings_.dirichletBcRight[0];
                //cdiscretization_->v(cdiscretization_->vJBegin(), j) = settings_.dirichletBcLeft[1];
                //cdiscretization_->v(cdiscretization_->vJEnd() - 1, j) = settings_.dirichletBcRight[1];
            }
            for (int i = cdiscretization_->pIBegin(); i < cdiscretization_->pIEnd(); i++)
            {
                for (int j = cdiscretization_->pJBegin(); j < cdiscretization_->pJEnd(); j++)
                {
                    for (int k = cdiscretization_->pKBegin();k<cdiscretization_->pKEnd();k++){
                        cdiscretization_->rho(i, j,k) = 1.;
                    }
                }
            }            

        }

        Collision();
        Streaming();
        // FillVelocitiesAndPressure();
        // FillVelocitiesAndPressure();
        LBMapplyBoundaryValues();

        FillVelocitiesAndPressure();
        // std::cout<<"*******************"<<std::endl;

        // LBMBounceBack();

        // std::cout<<"-------------------"<<std::endl;
        // FillVelocitiesAndPressure();
        // std::cout<<"+++++++++++++++++++"<<std::endl;
        outputWriterText_->writeFile(time);
        outputWriterParaview_->writeFile(time);

        // std::cout<<"                 "<<std::endl;
        // for(int k = 0; k<9;k++){
        //     std::cout<<"  dt:  "<<dt_<<"  "<<"########pdf0/2: "<<cdiscretization_->pdf(0,4,k)<<" "<<"           ########pdf0/1: "<<cdiscretization_->pdf(1,2,k)<<" "<<"            ########pdf0/0: "<<cdiscretization_->pdf(1,3,k)<<std::endl;
        // }

        time += dt_;

        //}
    }
};

void Computation::LBMapplyBoundaryValues()
{
    // top condition
    int pdfJEnd = cdiscretization_->pdfJEnd() - 1;
    int pdfJBegin = cdiscretization_->pdfJBegin();
    int pdfIBegin = cdiscretization_->pdfIBegin();
    int pdfIEnd = cdiscretization_->pdfIEnd() - 1;
    int pdfKBegin = cdiscretization_->pdfKBegin();
    int pdfKEnd  = cdiscretization_->pdfKEnd() - 1;

    //set top inner
    for(int i = pdfIBegin+1; i<pdfIEnd;i++)
    {

            for(int j = pdfJBegin+1;j<pdfJEnd;j++)
            {
                double rho_N = 1./(1.+settings_.dirichletBcTop[2])*(cdiscretization_->pdf(i,j,pdfKEnd,0)+cdiscretization_->pdf(i,j,pdfKEnd,1)+cdiscretization_->pdf(i,j,pdfKEnd,2)
                                                                    +cdiscretization_->pdf(i,j,pdfKEnd,3)+cdiscretization_->pdf(i,j,pdfKEnd,4)+cdiscretization_->pdf(i,j,pdfKEnd,7)
                                                                    +cdiscretization_->pdf(i,j,pdfKEnd,8)+cdiscretization_->pdf(i,j,pdfKEnd,9)+cdiscretization_->pdf(i,j,pdfKEnd,10)+
                                                                    2.*(cdiscretization_->pdf(i,j,pdfKEnd,5)+cdiscretization_->pdf(i,j,pdfKEnd,11)+cdiscretization_->pdf(i,j,pdfKEnd,14)+
                                                                    cdiscretization_->pdf(i,j,pdfKEnd,15)+cdiscretization_->pdf(i,j,pdfKEnd,18)));
                double N_x = 1./2.*(cdiscretization_->pdf(i,j,pdfKEnd,1)+cdiscretization_->pdf(i,j,pdfKEnd,7)+cdiscretization_->pdf(i,j,pdfKEnd,10)-(cdiscretization_->pdf(i,j,pdfKEnd,2)
                                                            +cdiscretization_->pdf(i,j,pdfKEnd,8)+cdiscretization_->pdf(i,j,pdfKEnd,9)))-1./3.*rho_N*settings_.dirichletBcTop[0];                                                            
                double N_y = 1./2.*(cdiscretization_->pdf(i,j,pdfKEnd,3)+cdiscretization_->pdf(i,j,pdfKEnd,7)+cdiscretization_->pdf(i,j,pdfKEnd,8)-(cdiscretization_->pdf(i,j,pdfKEnd,4)
                                                        +cdiscretization_->pdf(i,j,pdfKEnd,10)+cdiscretization_->pdf(i,j,pdfKEnd,9)))-1./3.*rho_N*settings_.dirichletBcTop[1]; 
                cdiscretization_->pdf(i,j,pdfKEnd,6) = cdiscretization_->pdf(i,j,pdfKEnd,5)-1./3.*rho_N*settings_.dirichletBcTop[2];
                cdiscretization_->pdf(i,j,pdfKEnd,12) = cdiscretization_->pdf(i,j,pdfKEnd,14)+rho_N/6.*(-settings_.dirichletBcTop[2]+settings_.dirichletBcTop[0])-N_x;
                cdiscretization_->pdf(i,j,pdfKEnd,13) = cdiscretization_->pdf(i,j,pdfKEnd,11)+rho_N/6.*(-settings_.dirichletBcTop[2]-settings_.dirichletBcTop[0])+N_x;
                cdiscretization_->pdf(i,j,pdfKEnd,16) = cdiscretization_->pdf(i,j,pdfKEnd,18)+rho_N/6.*(-settings_.dirichletBcTop[2]+settings_.dirichletBcTop[1])-N_y;
                cdiscretization_->pdf(i,j,pdfKEnd,17) = cdiscretization_->pdf(i,j,pdfKEnd,15)+rho_N/6.*(-settings_.dirichletBcTop[2]-settings_.dirichletBcTop[1])+N_y;

            }
        
    }

    //set bot inner
    for(int i = pdfIBegin+1; i<pdfIEnd;i++)
    {

            for(int j = pdfJBegin+1;j<pdfJEnd;j++)
            {
                double rho_N = 1./(1.-settings_.dirichletBcBottom[2])*(cdiscretization_->pdf(i,j,pdfKBegin,0)+cdiscretization_->pdf(i,j,pdfKBegin,1)+cdiscretization_->pdf(i,j,pdfKBegin,2)
                                                                    +cdiscretization_->pdf(i,j,pdfKBegin,3)+cdiscretization_->pdf(i,j,pdfKBegin,4)+cdiscretization_->pdf(i,j,pdfKBegin,7)
                                                                    +cdiscretization_->pdf(i,j,pdfKBegin,8)+cdiscretization_->pdf(i,j,pdfKBegin,9)+cdiscretization_->pdf(i,j,pdfKBegin,10)+
                                                                    2.*(cdiscretization_->pdf(i,j,pdfKBegin,6)+cdiscretization_->pdf(i,j,pdfKBegin,12)+cdiscretization_->pdf(i,j,pdfKBegin,13)+
                                                                    cdiscretization_->pdf(i,j,pdfKBegin,16)+cdiscretization_->pdf(i,j,pdfKBegin,17)));
                double N_x = 1./2.*(cdiscretization_->pdf(i,j,pdfKBegin,1)+cdiscretization_->pdf(i,j,pdfKBegin,7)+cdiscretization_->pdf(i,j,pdfKBegin,10)-(cdiscretization_->pdf(i,j,pdfKBegin,2)
                                                            +cdiscretization_->pdf(i,j,pdfKBegin,8)+cdiscretization_->pdf(i,j,pdfKBegin,9)))-1./3.*rho_N*settings_.dirichletBcBottom[0];                                                            
                double N_y = 1./2.*(cdiscretization_->pdf(i,j,pdfKBegin,3)+cdiscretization_->pdf(i,j,pdfKBegin,7)+cdiscretization_->pdf(i,j,pdfKBegin,8)-(cdiscretization_->pdf(i,j,pdfKBegin,4)
                                                        +cdiscretization_->pdf(i,j,pdfKBegin,10)+cdiscretization_->pdf(i,j,pdfKBegin,9)))-1./3.*rho_N*settings_.dirichletBcBottom[1];  
                cdiscretization_->pdf(i,j,pdfKBegin,5) = cdiscretization_->pdf(i,j,pdfKBegin,6)+1./3.*rho_N*settings_.dirichletBcBottom[2];
                cdiscretization_->pdf(i,j,pdfKBegin,11) = cdiscretization_->pdf(i,j,pdfKBegin,13)+rho_N/6.*(settings_.dirichletBcBottom[2]+settings_.dirichletBcBottom[0])-N_x;
                cdiscretization_->pdf(i,j,pdfKBegin,14) = cdiscretization_->pdf(i,j,pdfKBegin,12)+rho_N/6.*(settings_.dirichletBcBottom[2]-settings_.dirichletBcBottom[0])+N_x;
                cdiscretization_->pdf(i,j,pdfKBegin,15) = cdiscretization_->pdf(i,j,pdfKBegin,17)+rho_N/6.*(settings_.dirichletBcBottom[2]+settings_.dirichletBcBottom[1])-N_y;
                cdiscretization_->pdf(i,j,pdfKBegin,18) = cdiscretization_->pdf(i,j,pdfKBegin,16)+rho_N/6.*(settings_.dirichletBcBottom[2]-settings_.dirichletBcBottom[1])+N_y;

            }
        
    }

    //set left inner
    for(int j = pdfJBegin+1;j<pdfJEnd;j++)
    {
        for(int k = pdfKBegin+1;k<pdfKEnd;k++)
        {
            double rho_N = 1./(1.+settings_.dirichletBcLeft[0])*(cdiscretization_->pdf(pdfIBegin,j,k,0)+cdiscretization_->pdf(pdfIBegin,j,k,3)+cdiscretization_->pdf(pdfIBegin,j,k,15)
                                                                +cdiscretization_->pdf(pdfIBegin,j,k,5)+cdiscretization_->pdf(pdfIBegin,j,k,18)+cdiscretization_->pdf(pdfIBegin,j,k,4)
                                                                +cdiscretization_->pdf(pdfIBegin,j,k,17)+cdiscretization_->pdf(pdfIBegin,j,k,6)+cdiscretization_->pdf(pdfIBegin,j,k,16)+
                                                                2.*(cdiscretization_->pdf(pdfIBegin,j,k,2)+cdiscretization_->pdf(pdfIBegin,j,k,9)+cdiscretization_->pdf(pdfIBegin,j,k,8)
                                                                    +cdiscretization_->pdf(pdfIBegin,j,k,14)+cdiscretization_->pdf(pdfIBegin,j,k,13)));
            double N_y = 1./2.*(cdiscretization_->pdf(pdfIBegin,j,k,3)+cdiscretization_->pdf(pdfIBegin,j,k,15)+cdiscretization_->pdf(pdfIBegin,j,k,16)-(cdiscretization_->pdf(pdfIBegin,j,k,4)
                                        +cdiscretization_->pdf(pdfIBegin,j,k,18)+cdiscretization_->pdf(pdfIBegin,j,k,17)))-1./3.*rho_N*settings_.dirichletBcLeft[1];                                       
            double N_z = 1./2.*(cdiscretization_->pdf(pdfIBegin,j,k,5)+cdiscretization_->pdf(pdfIBegin,j,k,8)+cdiscretization_->pdf(pdfIBegin,j,k,15)-
                                                    (cdiscretization_->pdf(pdfIBegin,j,k,6)+cdiscretization_->pdf(pdfIBegin,j,k,16)+cdiscretization_->pdf(pdfIBegin,j,k,17)))-1./3.*rho_N*settings_.dirichletBcLeft[2];
            cdiscretization_->pdf(pdfIBegin,j,k,1) = cdiscretization_->pdf(pdfIBegin,j,k,2)+1./3.*rho_N*settings_.dirichletBcLeft[0];
            cdiscretization_->pdf(pdfIBegin,j,k,10) = cdiscretization_->pdf(pdfIBegin,j,k,8)+rho_N/6.*(settings_.dirichletBcLeft[0]-settings_.dirichletBcLeft[1])+N_y;
            cdiscretization_->pdf(pdfIBegin,j,k,7) = cdiscretization_->pdf(pdfIBegin,j,k,9)+rho_N/6.*(settings_.dirichletBcLeft[0]+settings_.dirichletBcLeft[1])-N_y;  
            cdiscretization_->pdf(pdfIBegin,j,k,11) = cdiscretization_->pdf(pdfIBegin,j,k,13)+rho_N/6.*(settings_.dirichletBcLeft[0]+settings_.dirichletBcLeft[2])-N_z;
            cdiscretization_->pdf(pdfIBegin,j,k,12) = cdiscretization_->pdf(pdfIBegin,j,k,14)+rho_N/6.*(settings_.dirichletBcLeft[0]-settings_.dirichletBcLeft[2])+N_z;            

        }
    }

    //set right inner
    for(int j = pdfJBegin+1;j<pdfJEnd;j++)
    {
        for(int k = pdfKBegin+1;k<pdfKEnd;k++)
        {
            double rho_N = 1./(1.-settings_.dirichletBcRight[0])*(cdiscretization_->pdf(pdfIBegin,j,k,0)+cdiscretization_->pdf(pdfIBegin,j,k,3)+cdiscretization_->pdf(pdfIBegin,j,k,15)
                                                                +cdiscretization_->pdf(pdfIBegin,j,k,5)+cdiscretization_->pdf(pdfIBegin,j,k,18)+cdiscretization_->pdf(pdfIBegin,j,k,4)
                                                                +cdiscretization_->pdf(pdfIBegin,j,k,17)+cdiscretization_->pdf(pdfIBegin,j,k,6)+cdiscretization_->pdf(pdfIBegin,j,k,16)+
                                                                2.*(cdiscretization_->pdf(pdfIBegin,j,k,1)+cdiscretization_->pdf(pdfIBegin,j,k,7)+cdiscretization_->pdf(pdfIBegin,j,k,11)
                                                                    +cdiscretization_->pdf(pdfIBegin,j,k,10)+cdiscretization_->pdf(pdfIBegin,j,k,12)));
            double N_y = 1./2.*(cdiscretization_->pdf(pdfIBegin,j,k,3)+cdiscretization_->pdf(pdfIBegin,j,k,15)+cdiscretization_->pdf(pdfIBegin,j,k,16)-(cdiscretization_->pdf(pdfIBegin,j,k,4)
                                        +cdiscretization_->pdf(pdfIBegin,j,k,18)+cdiscretization_->pdf(pdfIBegin,j,k,17)))-1./3.*rho_N*settings_.dirichletBcRight[1];                                       
            double N_z = 1./2.*(cdiscretization_->pdf(pdfIBegin,j,k,5)+cdiscretization_->pdf(pdfIBegin,j,k,8)+cdiscretization_->pdf(pdfIBegin,j,k,15)-
                                                    (cdiscretization_->pdf(pdfIBegin,j,k,6)+cdiscretization_->pdf(pdfIBegin,j,k,16)+cdiscretization_->pdf(pdfIBegin,j,k,17)))-1./3.*rho_N*settings_.dirichletBcRight[2];
            cdiscretization_->pdf(pdfIBegin,j,k,2) = cdiscretization_->pdf(pdfIBegin,j,k,1)-1./3.*rho_N*settings_.dirichletBcRight[0];
            cdiscretization_->pdf(pdfIBegin,j,k,8) = cdiscretization_->pdf(pdfIBegin,j,k,10)-rho_N/6.*(settings_.dirichletBcRight[0]-settings_.dirichletBcRight[1])-N_y;
            cdiscretization_->pdf(pdfIBegin,j,k,9) = cdiscretization_->pdf(pdfIBegin,j,k,7)-rho_N/6.*(settings_.dirichletBcRight[0]+settings_.dirichletBcRight[1])+N_y;  
            cdiscretization_->pdf(pdfIBegin,j,k,13) = cdiscretization_->pdf(pdfIBegin,j,k,11)-rho_N/6.*(settings_.dirichletBcRight[0]+settings_.dirichletBcRight[2])+N_z;
            cdiscretization_->pdf(pdfIBegin,j,k,14) = cdiscretization_->pdf(pdfIBegin,j,k,12)-rho_N/6.*(settings_.dirichletBcRight[0]-settings_.dirichletBcRight[2])-N_z;            

        }
    }

    //set front inner
    for(int i = pdfIBegin+1;i<pdfIEnd;i++)
    {
        for(int k = pdfKBegin+1;k<pdfKEnd;k++)
        {
            double rho_N = 1./(1.-settings_.dirichletBcFront[1])*(cdiscretization_->pdf(i,pdfJBegin,k,0)+cdiscretization_->pdf(i,pdfJBegin,k,1)+cdiscretization_->pdf(i,pdfJBegin,k,2)
                                                                +cdiscretization_->pdf(i,pdfJBegin,k,5)+cdiscretization_->pdf(i,pdfJBegin,k,6)+cdiscretization_->pdf(i,pdfJBegin,k,12)
                                                                +cdiscretization_->pdf(i,pdfJBegin,k,11)+cdiscretization_->pdf(i,pdfJBegin,k,14)+cdiscretization_->pdf(i,pdfJBegin,k,13)+
                                                                2.*(cdiscretization_->pdf(i,pdfJBegin,k,4)+cdiscretization_->pdf(i,pdfJBegin,k,9)+cdiscretization_->pdf(i,pdfJBegin,k,18)
                                                                    +cdiscretization_->pdf(i,pdfJBegin,k,10)+cdiscretization_->pdf(i,pdfJBegin,k,17)));
            double N_x = 1./2.*(cdiscretization_->pdf(i,pdfJBegin,k,1)+cdiscretization_->pdf(i,pdfJBegin,k,11)+cdiscretization_->pdf(i,pdfJBegin,k,12)-(cdiscretization_->pdf(i,pdfJBegin,k,2)
                                        +cdiscretization_->pdf(i,pdfJBegin,k,14)+cdiscretization_->pdf(i,pdfJBegin,k,13)))-1./3.*rho_N*settings_.dirichletBcFront[0];                                       
            double N_y = 1./2.*(cdiscretization_->pdf(i,pdfJBegin,k,3)+cdiscretization_->pdf(i,pdfJBegin,k,15)+cdiscretization_->pdf(i,pdfJBegin,k,16)-(cdiscretization_->pdf(i,pdfJBegin,k,4)
                                        +cdiscretization_->pdf(i,pdfJBegin,k,18)+cdiscretization_->pdf(i,pdfJBegin,k,17)))-1./3.*rho_N*settings_.dirichletBcFront[1];
            double N_z = 1./2.*(cdiscretization_->pdf(i,pdfJBegin,k,5)+cdiscretization_->pdf(i,pdfJBegin,k,11)+cdiscretization_->pdf(i,pdfJBegin,k,14)-
                                                    (cdiscretization_->pdf(i,pdfJBegin,k,6)+cdiscretization_->pdf(i,pdfJBegin,k,12)+cdiscretization_->pdf(i,pdfJBegin,k,13)))
                                                    -1./3.*rho_N*settings_.dirichletBcFront[2];
            cdiscretization_->pdf(i,pdfJBegin,k,3) = cdiscretization_->pdf(i,pdfJBegin,k,4)+1./3.*rho_N*settings_.dirichletBcFront[1];
            cdiscretization_->pdf(i,pdfJBegin,k,7) = cdiscretization_->pdf(i,pdfJBegin,k,9)+rho_N/6.*(settings_.dirichletBcFront[0]+settings_.dirichletBcFront[1])-N_x;
            cdiscretization_->pdf(i,pdfJBegin,k,8) = cdiscretization_->pdf(i,pdfJBegin,k,10)+rho_N/6.*(-settings_.dirichletBcFront[0]+settings_.dirichletBcFront[1])+N_x;  
            cdiscretization_->pdf(i,pdfJBegin,k,15) = cdiscretization_->pdf(i,pdfJBegin,k,17)+rho_N/6.*(settings_.dirichletBcFront[1]+settings_.dirichletBcFront[2])-N_z;
            cdiscretization_->pdf(i,pdfJBegin,k,16) = cdiscretization_->pdf(i,pdfJBegin,k,18)+rho_N/6.*(settings_.dirichletBcFront[1]-settings_.dirichletBcFront[2])+N_y;            

        }
    }

    //set back inner
    for(int i = pdfIBegin+1;i<pdfIEnd;i++)
    {
        for(int k = pdfKBegin+1;k<pdfKEnd;k++)
        {
            double rho_N = 1./(1.+settings_.dirichletBcBack[1])*(cdiscretization_->pdf(i,pdfJEnd,k,0)+cdiscretization_->pdf(i,pdfJEnd,k,1)+cdiscretization_->pdf(i,pdfJEnd,k,2)
                                                                +cdiscretization_->pdf(i,pdfJEnd,k,5)+cdiscretization_->pdf(i,pdfJEnd,k,6)+cdiscretization_->pdf(i,pdfJEnd,k,12)
                                                                +cdiscretization_->pdf(i,pdfJEnd,k,11)+cdiscretization_->pdf(i,pdfJEnd,k,14)+cdiscretization_->pdf(i,pdfJEnd,k,13)+
                                                                2.*(cdiscretization_->pdf(i,pdfJEnd,k,3)+cdiscretization_->pdf(i,pdfJEnd,k,7)+cdiscretization_->pdf(i,pdfJEnd,k,8)
                                                                    +cdiscretization_->pdf(i,pdfJEnd,k,16)+cdiscretization_->pdf(i,pdfJEnd,k,15)));
            double N_x = 1./2.*(cdiscretization_->pdf(i,pdfJEnd,k,1)+cdiscretization_->pdf(i,pdfJEnd,k,11)+cdiscretization_->pdf(i,pdfJEnd,k,12)-(cdiscretization_->pdf(i,pdfJEnd,k,2)
                                        +cdiscretization_->pdf(i,pdfJEnd,k,14)+cdiscretization_->pdf(i,pdfJEnd,k,13)))-1./3.*rho_N*settings_.dirichletBcBack[0];                                       
            double N_y = 1./2.*(cdiscretization_->pdf(i,pdfJBegin,k,3)+cdiscretization_->pdf(i,pdfJBegin,k,15)+cdiscretization_->pdf(i,pdfJBegin,k,16)-(cdiscretization_->pdf(i,pdfJBegin,k,4)
                                        +cdiscretization_->pdf(i,pdfJBegin,k,18)+cdiscretization_->pdf(i,pdfJBegin,k,17)))-1./3.*rho_N*settings_.dirichletBcBack[1];
            double N_z = 1./2.*(cdiscretization_->pdf(i,pdfJEnd,k,5)+cdiscretization_->pdf(i,pdfJEnd,k,11)+cdiscretization_->pdf(i,pdfJEnd,k,14)-
                                                    (cdiscretization_->pdf(i,pdfJEnd,k,6)+cdiscretization_->pdf(i,pdfJEnd,k,12)+cdiscretization_->pdf(i,pdfJEnd,k,13)))
                                                    -1./3.*rho_N*settings_.dirichletBcBack[2];
            cdiscretization_->pdf(i,pdfJEnd,k,4) = cdiscretization_->pdf(i,pdfJEnd,k,3)-1./3.*rho_N*settings_.dirichletBcBack[1];
            cdiscretization_->pdf(i,pdfJEnd,k,9) = cdiscretization_->pdf(i,pdfJEnd,k,7)-rho_N/6.*(settings_.dirichletBcBack[0]+settings_.dirichletBcBack[1])+N_x;
            cdiscretization_->pdf(i,pdfJEnd,k,10) = cdiscretization_->pdf(i,pdfJEnd,k,8)-rho_N/6.*(-settings_.dirichletBcBack[0]+settings_.dirichletBcBack[1])-N_x;  
            cdiscretization_->pdf(i,pdfJEnd,k,17) = cdiscretization_->pdf(i,pdfJEnd,k,15)-rho_N/6.*(settings_.dirichletBcBack[1]+settings_.dirichletBcBack[2])+N_z;
            cdiscretization_->pdf(i,pdfJEnd,k,18) = cdiscretization_->pdf(i,pdfJEnd,k,16)-rho_N/6.*(settings_.dirichletBcBack[1]-settings_.dirichletBcBack[2])-N_y;            

        }
    }

// //top left
    // cdiscretization_->u(0,pdfJEnd) = 0.0;
    // cdiscretization_->v(0,pdfJEnd) = 0.0;
    // cdiscretization_->pdf(0,pdfJEnd,8) = cdiscretization_->pdf(0,pdfJEnd,6);
    // cdiscretization_->pdf(0,pdfJEnd,4) = cdiscretization_->pdf(0,pdfJEnd,2);
    // cdiscretization_->pdf(0,pdfJEnd,1) = cdiscretization_->pdf(0,pdfJEnd,3);
    // cdiscretization_->pdf(0,pdfJEnd,7) = 0.5*(cdiscretization_->rho(0,pdfJEnd)-(cdiscretization_->pdf(0,pdfJEnd,0)+cdiscretization_->pdf(0,pdfJEnd,1)+cdiscretization_->pdf(0,pdfJEnd,2)
    //                                         +cdiscretization_->pdf(0,pdfJEnd,3)+cdiscretization_->pdf(0,pdfJEnd,4)+cdiscretization_->pdf(0,pdfJEnd,6)+cdiscretization_->pdf(0,pdfJEnd,8)));
    // cdiscretization_->pdf(0,pdfJEnd,5) = cdiscretization_->pdf(0,pdfJEnd,7);
    // //top right
    // cdiscretization_->u(pdfIEnd,pdfJEnd) = 0.0;
    // cdiscretization_->v(pdfIEnd,pdfJEnd) = 0.0;
    // cdiscretization_->pdf(pdfIEnd,pdfJEnd,7) = cdiscretization_->pdf(pdfIEnd,pdfJEnd,5);
    // cdiscretization_->pdf(pdfIEnd,pdfJEnd,4) = cdiscretization_->pdf(pdfIEnd,pdfJEnd,2);
    // cdiscretization_->pdf(pdfIEnd,pdfJEnd,1) = cdiscretization_->pdf(pdfIEnd,pdfJEnd,3);
    // cdiscretization_->pdf(pdfIEnd,pdfJEnd,8) = 0.5*(cdiscretization_->rho(pdfIEnd,pdfJEnd)-(cdiscretization_->pdf(pdfIEnd,pdfJEnd,0)
    //                                         +cdiscretization_->pdf(pdfIEnd,pdfJEnd,1)+cdiscretization_->pdf(pdfIEnd,pdfJEnd,2)
    //                                         +cdiscretization_->pdf(pdfIEnd,pdfJEnd,3)+cdiscretization_->pdf(pdfIEnd,pdfJEnd,4)
    //                                         +cdiscretization_->pdf(pdfIEnd,pdfJEnd,5)+cdiscretization_->pdf(pdfIEnd,pdfJEnd,7)));
    // cdiscretization_->pdf(pdfIEnd,pdfJEnd,6) = cdiscretization_->pdf(pdfIEnd,pdfJEnd,8);

    // //bottom right
    // cdiscretization_->u(pdfIEnd,0) = 0.0;
    // cdiscretization_->v(pdfIEnd,0) = 0.0;
    // cdiscretization_->pdf(pdfIEnd,0,6) = cdiscretization_->pdf(pdfIEnd,0,8);
    // cdiscretization_->pdf(pdfIEnd,0,2) = cdiscretization_->pdf(pdfIEnd,0,4);
    // cdiscretization_->pdf(pdfIEnd,0,3) = cdiscretization_->pdf(pdfIEnd,0,1);
    // cdiscretization_->pdf(pdfIEnd,0,7) = 0.5*(cdiscretization_->rho(pdfIEnd,0)-(cdiscretization_->pdf(pdfIEnd,0,0)+cdiscretization_->pdf(pdfIEnd,0,1)+cdiscretization_->pdf(pdfIEnd,0,2)
    //                                         +cdiscretization_->pdf(pdfIEnd,0,3)+cdiscretization_->pdf(pdfIEnd,0,4)+cdiscretization_->pdf(pdfIEnd,0,8)+cdiscretization_->pdf(pdfIEnd,0,6)));
    // cdiscretization_->pdf(pdfIEnd,0,5) = cdiscretization_->pdf(pdfIEnd,0,7);

    // //bottom left
    // cdiscretization_->u(0,0) = 0.0;
    // cdiscretization_->v(0,0) = 0.0;
    // cdiscretization_->pdf(0,0,5) = cdiscretization_->pdf(0,0,7);
    // cdiscretization_->pdf(0,0,2) = cdiscretization_->pdf(0,0,4);
    // cdiscretization_->pdf(0,0,1) = cdiscretization_->pdf(0,0,3);
    // cdiscretization_->pdf(0,0,8) = 0.5*(cdiscretization_->rho(0,0)-(cdiscretization_->pdf(0,0,0)+cdiscretization_->pdf(0,0,1)+cdiscretization_->pdf(0,0,2)
    //                                         +cdiscretization_->pdf(0,0,3)+cdiscretization_->pdf(0,0,4)+cdiscretization_->pdf(0,0,5)+cdiscretization_->pdf(0,0,7)));
    // cdiscretization_->pdf(0,0,6) = cdiscretization_->pdf(0,0,8);
};

void Computation::LBMBounceBack(){
    //    top condition

};

void Computation::FillVelocitiesAndPressure()
{
    for (int i = cdiscretization_->pIBegin(); i < cdiscretization_->pIEnd(); i++)
    {
        for (int j = cdiscretization_->pJBegin(); j < cdiscretization_->pJEnd(); j++)
        {
            double p = 0.;
            double u = 0.;
            double v = 0.;
            double w = 0.;
            double rho = 0.;
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd(); k++)
            {
                for (int l = 0; l < 9; l++)
                {
                    rho += cdiscretization_->pdf(i, j, k, l);
                    p += cdiscretization_->pdf(i, j, k, l) * cs_ * cs_;
                    u += cdiscretization_->pdf(i, j, k, l) * cdiscretization_->ci_x(l);
                    v += cdiscretization_->pdf(i, j, k, l) * cdiscretization_->ci_y(l);
                    w += cdiscretization_->pdf(i, j, k, l) * cdiscretization_->ci_z(l);
                    // std::cout<<"###### "<<cdiscretization_->pdf(i,j,k)<< "####rho "<<cdiscretization_->rho(i,j)<<std::endl;
                }
            

            cdiscretization_->rho(i, j, k) = (double)rho;
            cdiscretization_->p(i, j, k) = (double)p;
            cdiscretization_->u(i, j, k) = (double)u / rho;
            cdiscretization_->v(i, j, k) = (double)v / rho;
            cdiscretization_->w(i,j,k) = (double) w / rho;
            }
        }
    }
};

void Computation::Collision()
{
    int pdfJEnd = cdiscretization_->pdfJEnd();
    int pdfJBegin = cdiscretization_->pdfJBegin();
    int pdfIBegin = cdiscretization_->pdfIBegin();
    int pdfIEnd = cdiscretization_->pdfIEnd();

    for (int i = pdfIBegin; i < pdfIEnd; i++)
    {
        for (int j = pdfJBegin; j < pdfJEnd; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd(); k++)
            {
                for (int l = 0; l < 19; l++)
                {
                    cdiscretization_->pdfeq(i, j, k, l) = cdiscretization_->wi(l) * cdiscretization_->rho(i, j, k) * (1. + 
                    (cdiscretization_->u(i, j, k) * cdiscretization_->ci_x(l) + cdiscretization_->v(i, j, k) * cdiscretization_->ci_y(l) + cdiscretization_->w(i, j, k) * cdiscretization_->ci_z(l)) / (cs_ * cs_) +
                     1. / 2. * ((cdiscretization_->u(i, j, k) * cdiscretization_->ci_x(l) + cdiscretization_->v(i, j, k) * cdiscretization_->ci_y(l) + cdiscretization_->w(i, j, k) * cdiscretization_->ci_z(l)) / (cs_ * cs_)) 
                     * ((cdiscretization_->u(i, j, k) * cdiscretization_->ci_x(l) + cdiscretization_->v(i, j, k) * cdiscretization_->ci_y(l) + cdiscretization_->w(i, j, k) * cdiscretization_->ci_z(l)) / (cs_ * cs_)) 
                     - 1. / 2. * (cdiscretization_->u(i, j, k) * cdiscretization_->u(i, j, k) + cdiscretization_->v(i, j, k) * cdiscretization_->v(i, j, k) + cdiscretization_->w(i, j, k) * cdiscretization_->w(i, j, k)) / (cs_ * cs_));
                    // std::cout<<"###### "<<cdiscretization_->pdfeq(i,j,k)<< "####rho "<<cdiscretization_->rho(i,j)<<std::endl;
                }
            }
        }
    }
};

void Computation::Streaming()
{
    // write pdf in pdfold
    for (int i = cdiscretization_->pdfIBegin(); i < cdiscretization_->pdfIEnd(); i++)
    {
        for (int j = cdiscretization_->pdfJBegin(); j < cdiscretization_->pdfJEnd(); j++)
        {
            for (int l = cdiscretization_->pdfKBegin(); l < cdiscretization_->pdfKEnd(); l++)
            {
                for (int k = 0; k < 19; k++)
                {
                    cdiscretization_->pdf(i, j, l, k) = cdiscretization_->pdf(i, j, l, k) - dt_ / tau_ * (cdiscretization_->pdf(i, j, l, k) - cdiscretization_->pdfeq(i, j, l, k));
                    cdiscretization_->pdfold(i, j, l, k) = cdiscretization_->pdf(i, j, l, k);
                }
            }
        }
    }
    int pdfIBeginInner = cdiscretization_->pdfIBegin() + 1;
    int pdfIEndInner = cdiscretization_->pdfIEnd() - 1;
    int pdfJBeginInner = cdiscretization_->pdfJBegin() + 1;
    int pdfJEndInner = cdiscretization_->pdfJEnd() - 1;

    for (int i = pdfIBeginInner - 1; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd(); k++)
            {
                cdiscretization_->pdf(i, j, k, 0) = cdiscretization_->pdfold(i, j, k, 0);
            }
        }
    }
    for (int i = pdfIBeginInner - 1; i < pdfIEndInner; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd(); k++){
                cdiscretization_->pdf(i + 1, j,k, 1) = cdiscretization_->pdfold(i, j,k, 1);
            }
        }
    }
    for (int i = cdiscretization_->pdfIBegin() + 1; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd(); k++){
                cdiscretization_->pdf(i - 1, j, k, 2) = cdiscretization_->pdfold(i, j, k, 2);
            }
        }
    }

    for (int i = pdfIBeginInner - 1; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd(); k++)
            {
                cdiscretization_->pdf(i, j + 1, k, 3) = cdiscretization_->pdfold(i, j, k, 3);
            }
        }
    }
    for (int i = pdfIBeginInner - 1; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd(); k++)
            {
                cdiscretization_->pdf(i, j - 1,k, 4) = cdiscretization_->pdfold(i, j,k, 4);
            }
        }
    }

    for (int i = pdfIBeginInner - 1; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd() - 1; k++)
            {
                cdiscretization_->pdf(i, j, k + 1, 5) = cdiscretization_->pdfold(i, j, k, 5);
            }
        }
    }

    for (int i = pdfIBeginInner - 1; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin() + 1; k < cdiscretization_->pdfKEnd(); k++)
            {
                cdiscretization_->pdf(i, j, k - 1, 6) = cdiscretization_->pdfold(i, j, k, 6);
            }
        }
    }
    for (int i = pdfIBeginInner - 1; i < pdfIEndInner; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd(); k++)
            {
                cdiscretization_->pdf(i + 1, j + 1, k, 7) = cdiscretization_->pdfold(i, j, k, 7);
            }
        }
    }
    for (int i = pdfIBeginInner; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd(); k++)
            {
                cdiscretization_->pdf(i - 1, j + 1, k, 8) = cdiscretization_->pdfold(i, j, k, 8);
            }
        }
    }

    for (int i = pdfIBeginInner; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd(); k++)
            {
                cdiscretization_->pdf(i - 1, j - 1, k, 9) = cdiscretization_->pdfold(i, j, k, 9);
            }
        }
    }

    for (int i = pdfIBeginInner + 1; i < pdfIEndInner; i++)
    {
        for (int j = pdfJBeginInner; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd(); k++)
            {
                cdiscretization_->pdf(i + 1, j - 1, k, 10) = cdiscretization_->pdfold(i, j, k, 10);
            }
        }
    }

    for (int i = pdfIBeginInner - 1; i < pdfIEndInner; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd() - 1; k++)
            {
                cdiscretization_->pdf(i + 1, j, k + 1, 11) = cdiscretization_->pdfold(i, j, k, 11);
            }
        }
    }

    for (int i = pdfIBeginInner - 1; i < pdfIEndInner; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin() + 1; k < cdiscretization_->pdfKEnd(); k++)
            {
                cdiscretization_->pdf(i + 1, j, k - 1, 12) = cdiscretization_->pdfold(i, j, k, 12);
            }
        }
    }

    for (int i = pdfIBeginInner; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin() + 1; k < cdiscretization_->pdfKEnd(); k++)
            {
                cdiscretization_->pdf(i - 1, j, k - 1, 13) = cdiscretization_->pdfold(i, j, k, 13);
            }
        }
    }

    for (int i = pdfIBeginInner; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd() - 1; k++)
            {
                cdiscretization_->pdf(i - 1, j, k + 1, 14) = cdiscretization_->pdfold(i, j, k, 14);
            }
        }
    }

    for (int i = pdfIBeginInner - 1; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd() - 1; k++)
            {
                cdiscretization_->pdf(i, j + 1, k + 1, 15) = cdiscretization_->pdfold(i, j, k, 15);
            }
        }
    }

    for (int i = pdfIBeginInner - 1; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner - 1; j < pdfJEndInner; j++)
        {
            for (int k = cdiscretization_->pdfKBegin() + 1; k < cdiscretization_->pdfKEnd(); k++)
            {
                cdiscretization_->pdf(i, j + 1, k - 1, 16) = cdiscretization_->pdfold(i, j, k, 16);
            }
        }
    }

    for (int i = pdfIBeginInner - 1; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin() + 1; k < cdiscretization_->pdfKEnd(); k++)
            {
                cdiscretization_->pdf(i, j - 1, k - 1, 17) = cdiscretization_->pdfold(i, j, k, 17);
            }
        }
    }

    for (int i = pdfIBeginInner - 1; i < pdfIEndInner + 1; i++)
    {
        for (int j = pdfJBeginInner; j < pdfJEndInner + 1; j++)
        {
            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd() - 1; k++)
            {
                cdiscretization_->pdf(i, j - 1, k + 1, 18) = cdiscretization_->pdfold(i, j, k, 18);
            }
        }
    }
};
