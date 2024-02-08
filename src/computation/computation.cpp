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

    for (int i = 0; i < 3; i++)
        meshWidth_[i] = settings_.physicalSize[i] / settings_.nCells[i];

    // doLBM = true;
    cdiscretization_ = std::make_shared<LbmDiscretization>(settings_.nCells, meshWidth_);

    cs_ = (double)1. / sqrt(3);
    if(settings_.rightBcPressure){
        nu_ = settings_.dirichletBcLeft[0] * settings_.L_lbm / settings_.re;
    }else{
        nu_ = settings_.dirichletBcTop[0] * settings_.L_lbm / settings_.re;
    }
    tau_ = .5 + nu_ / (cs_ * cs_);
    dt_ = 1.; // 
    std::cout << dt_ << std::endl;


    std::cout << tau_ << std::endl;
    outputWriterText_ = std::make_unique<OutputWriterText>(cdiscretization_);
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(cdiscretization_);

};

/**
 * Run the whole simulation until tend
 */
void Computation::runSimulation()
{
    int t_iter = 0;
    double time = 0.0;
    double dtp_ = settings_.re*nu_/(settings_.nCells[0]*settings_.nCells[0]);
    while (time < settings_.endTime)
    {

        // settings_.dirichletBcTop[0] = 1000.*(1.0-exp(-t_iter*t_iter/(2.0*10*settings_.nCells[0]*10*settings_.nCells[0])));
        t_iter++;

        // std::cout<<"+++++++++++++++++++"<<std::endl;
        if (time == 0.0)
        {
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
      

        LBMapplyBoundaryValues();
        FillVelocitiesAndPressure();
        if(t_iter%settings_.deltawrite_==0||time ==0.0)
            outputWriterParaview_->writeFile(time); 

        time += dtp_;

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
                cdiscretization_->pdf(i,j,pdfKEnd,6) = cdiscretization_->pdf(i,j,pdfKEnd,5)-1./3.*rho_N*settings_.dirichletBcTop[2];
                cdiscretization_->pdf(i,j,pdfKEnd,12) = cdiscretization_->pdf(i,j,pdfKEnd,14)-1./4.*(cdiscretization_->pdf(i,j,pdfKEnd,1)-cdiscretization_->pdf(i,j,pdfKEnd,2))+rho_N/6.*(-settings_.dirichletBcTop[0]);
                cdiscretization_->pdf(i,j,pdfKEnd,13) = cdiscretization_->pdf(i,j,pdfKEnd,11)+1./4.*(cdiscretization_->pdf(i,j,pdfKEnd,1)-cdiscretization_->pdf(i,j,pdfKEnd,2))+rho_N/6.*(settings_.dirichletBcTop[0]);
                cdiscretization_->pdf(i,j,pdfKEnd,16) = cdiscretization_->pdf(i,j,pdfKEnd,18)-1./4.*(cdiscretization_->pdf(i,j,pdfKEnd,3)-cdiscretization_->pdf(i,j,pdfKEnd,4))+rho_N/6.*(-settings_.dirichletBcTop[2]);
                cdiscretization_->pdf(i,j,pdfKEnd,17) = cdiscretization_->pdf(i,j,pdfKEnd,15)+1./4.*(cdiscretization_->pdf(i,j,pdfKEnd,3)-cdiscretization_->pdf(i,j,pdfKEnd,4))+rho_N/6.*(-settings_.dirichletBcTop[2]);
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
                cdiscretization_->pdf(i,j,pdfKBegin,5) = cdiscretization_->pdf(i,j,pdfKBegin,6);
                cdiscretization_->pdf(i,j,pdfKBegin,11) = cdiscretization_->pdf(i,j,pdfKBegin,13);
                cdiscretization_->pdf(i,j,pdfKBegin,14) = cdiscretization_->pdf(i,j,pdfKBegin,12);
                cdiscretization_->pdf(i,j,pdfKBegin,15) = cdiscretization_->pdf(i,j,pdfKBegin,17);
                cdiscretization_->pdf(i,j,pdfKBegin,18) = cdiscretization_->pdf(i,j,pdfKBegin,16);

            }
        
    }

    //set left inner
    for(int j = pdfJBegin+1;j<pdfJEnd;j++)
    {
        for(int k = pdfKBegin+1;k<pdfKEnd;k++)
        {
            if(settings_.rightBcPressure){
            cdiscretization_->u(cdiscretization_->uIBegin(),j,k) = settings_.dirichletBcLeft[0];
            cdiscretization_->v(cdiscretization_->vIBegin(),j,k) = 0.0;
            cdiscretization_->w(cdiscretization_->wIBegin(),j,k) = 0.0;
            double rho_N = 1./(1.+settings_.dirichletBcLeft[0])*(cdiscretization_->pdf(pdfIBegin,j,k,0)+cdiscretization_->pdf(pdfIBegin,j,k,3)+cdiscretization_->pdf(pdfIBegin,j,k,15)
                                                                +cdiscretization_->pdf(pdfIBegin,j,k,5)+cdiscretization_->pdf(pdfIBegin,j,k,18)+cdiscretization_->pdf(pdfIBegin,j,k,4)
                                                                +cdiscretization_->pdf(pdfIBegin,j,k,17)+cdiscretization_->pdf(pdfIBegin,j,k,6)+cdiscretization_->pdf(pdfIBegin,j,k,16)+
                                                                2.*(cdiscretization_->pdf(pdfIBegin,j,k,2)+cdiscretization_->pdf(pdfIBegin,j,k,9)+cdiscretization_->pdf(pdfIBegin,j,k,8)
                                                                    +cdiscretization_->pdf(pdfIBegin,j,k,14)+cdiscretization_->pdf(pdfIBegin,j,k,13)));
            cdiscretization_->pdf(pdfIBegin,j,k,1) = cdiscretization_->pdf(pdfIBegin,j,k,2)+1./3.*rho_N*settings_.dirichletBcLeft[0];
            cdiscretization_->pdf(pdfIBegin,j,k,10) = cdiscretization_->pdf(pdfIBegin,j,k,8)+rho_N/6.*(settings_.dirichletBcLeft[0]);
            cdiscretization_->pdf(pdfIBegin,j,k,7) = cdiscretization_->pdf(pdfIBegin,j,k,9)+rho_N/6.*(settings_.dirichletBcLeft[0]);  
            cdiscretization_->pdf(pdfIBegin,j,k,11) = cdiscretization_->pdf(pdfIBegin,j,k,13)+rho_N/6.*(settings_.dirichletBcLeft[0]);
            cdiscretization_->pdf(pdfIBegin,j,k,12) = cdiscretization_->pdf(pdfIBegin,j,k,14)+rho_N/6.*(settings_.dirichletBcLeft[0]); 
            }else{
                double rho_N = 1./(1.+settings_.dirichletBcLeft[0])*(cdiscretization_->pdf(pdfIBegin,j,k,0)+cdiscretization_->pdf(pdfIBegin,j,k,3)+cdiscretization_->pdf(pdfIBegin,j,k,15)
                                                                    +cdiscretization_->pdf(pdfIBegin,j,k,5)+cdiscretization_->pdf(pdfIBegin,j,k,18)+cdiscretization_->pdf(pdfIBegin,j,k,4)
                                                                    +cdiscretization_->pdf(pdfIBegin,j,k,17)+cdiscretization_->pdf(pdfIBegin,j,k,6)+cdiscretization_->pdf(pdfIBegin,j,k,16)+
                                                                    2.*(cdiscretization_->pdf(pdfIBegin,j,k,2)+cdiscretization_->pdf(pdfIBegin,j,k,9)+cdiscretization_->pdf(pdfIBegin,j,k,8)
                                                                        +cdiscretization_->pdf(pdfIBegin,j,k,14)+cdiscretization_->pdf(pdfIBegin,j,k,13)));
                cdiscretization_->pdf(pdfIBegin,j,k,1) = cdiscretization_->pdf(pdfIBegin,j,k,2)+1./3.*rho_N*settings_.dirichletBcLeft[0];
                cdiscretization_->pdf(pdfIBegin,j,k,10) = cdiscretization_->pdf(pdfIBegin,j,k,8);
                cdiscretization_->pdf(pdfIBegin,j,k,7) = cdiscretization_->pdf(pdfIBegin,j,k,9);  
                cdiscretization_->pdf(pdfIBegin,j,k,11) = cdiscretization_->pdf(pdfIBegin,j,k,13);
                cdiscretization_->pdf(pdfIBegin,j,k,12) = cdiscretization_->pdf(pdfIBegin,j,k,14); 
            }
           

        }
    }

    //set right inner
    for(int j = pdfJBegin+1;j<pdfJEnd;j++)
    {
        for(int k = pdfKBegin+1;k<pdfKEnd;k++)
        {
            if(settings_.rightBcPressure){
                cdiscretization_->w(cdiscretization_->wIEnd()-1,j,k) = 0.0;
                cdiscretization_->v(cdiscretization_->vIEnd()-1,j,k) = 0.0;
                cdiscretization_->rho(cdiscretization_->rhoIEnd()-1,j,k) = settings_.rhoRight;
                cdiscretization_->u(cdiscretization_->uIEnd()-1,j,k) = 1.-(cdiscretization_->pdf(pdfIEnd,j,k,0)+cdiscretization_->pdf(pdfIEnd,j,k,3)+cdiscretization_->pdf(pdfIEnd,j,k,15)
                                                                    +cdiscretization_->pdf(pdfIEnd,j,k,5)+cdiscretization_->pdf(pdfIEnd,j,k,18)+cdiscretization_->pdf(pdfIEnd,j,k,4)
                                                                    +cdiscretization_->pdf(pdfIEnd,j,k,17)+cdiscretization_->pdf(pdfIEnd,j,k,6)+cdiscretization_->pdf(pdfIEnd,j,k,16)+
                                                                    2.*(cdiscretization_->pdf(pdfIEnd,j,k,1)+cdiscretization_->pdf(pdfIEnd,j,k,7)+cdiscretization_->pdf(pdfIEnd,j,k,11)
                                                                        +cdiscretization_->pdf(pdfIEnd,j,k,10)+cdiscretization_->pdf(pdfIEnd,j,k,12)))/settings_.rhoRight;
                cdiscretization_->pdf(pdfIEnd,j,k,2) = cdiscretization_->pdf(pdfIEnd,j,k,1)-1./3.*settings_.rhoRight*cdiscretization_->w(cdiscretization_->wIEnd()-1,j,k);
                cdiscretization_->pdf(pdfIEnd,j,k,8) = cdiscretization_->pdf(pdfIEnd,j,k,10)-cdiscretization_->rho(cdiscretization_->rhoIEnd()-1,j,k)/6.*cdiscretization_->w(cdiscretization_->wIEnd()-1,j,k);
                cdiscretization_->pdf(pdfIEnd,j,k,9) = cdiscretization_->pdf(pdfIEnd,j,k,7)-cdiscretization_->rho(cdiscretization_->rhoIEnd()-1,j,k)/6.*cdiscretization_->w(cdiscretization_->wIEnd()-1,j,k);  
                cdiscretization_->pdf(pdfIEnd,j,k,13) = cdiscretization_->pdf(pdfIEnd,j,k,11)-cdiscretization_->rho(cdiscretization_->rhoIEnd()-1,j,k)/6.*cdiscretization_->w(cdiscretization_->wIEnd()-1,j,k);
                cdiscretization_->pdf(pdfIEnd,j,k,14) = cdiscretization_->pdf(pdfIEnd,j,k,12)-cdiscretization_->rho(cdiscretization_->rhoIEnd()-1,j,k)/6.*cdiscretization_->w(cdiscretization_->wIEnd()-1,j,k); 
            }else{
                double rho_N = 1./(1.-settings_.dirichletBcRight[0])*(cdiscretization_->pdf(pdfIEnd,j,k,0)+cdiscretization_->pdf(pdfIEnd,j,k,3)+cdiscretization_->pdf(pdfIEnd,j,k,15)
                                                                    +cdiscretization_->pdf(pdfIEnd,j,k,5)+cdiscretization_->pdf(pdfIEnd,j,k,18)+cdiscretization_->pdf(pdfIEnd,j,k,4)
                                                                    +cdiscretization_->pdf(pdfIEnd,j,k,17)+cdiscretization_->pdf(pdfIEnd,j,k,6)+cdiscretization_->pdf(pdfIEnd,j,k,16)+
                                                                    2.*(cdiscretization_->pdf(pdfIEnd,j,k,1)+cdiscretization_->pdf(pdfIEnd,j,k,7)+cdiscretization_->pdf(pdfIEnd,j,k,11)
                                                                        +cdiscretization_->pdf(pdfIEnd,j,k,10)+cdiscretization_->pdf(pdfIEnd,j,k,12)));
                cdiscretization_->pdf(pdfIEnd,j,k,2) = cdiscretization_->pdf(pdfIEnd,j,k,1);
                cdiscretization_->pdf(pdfIEnd,j,k,8) = cdiscretization_->pdf(pdfIEnd,j,k,10);
                cdiscretization_->pdf(pdfIEnd,j,k,9) = cdiscretization_->pdf(pdfIEnd,j,k,7);
                cdiscretization_->pdf(pdfIEnd,j,k,13) = cdiscretization_->pdf(pdfIEnd,j,k,11);
                cdiscretization_->pdf(pdfIEnd,j,k,14) = cdiscretization_->pdf(pdfIEnd,j,k,12);
            }        

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
            cdiscretization_->pdf(i,pdfJBegin,k,3) = cdiscretization_->pdf(i,pdfJBegin,k,4);
            cdiscretization_->pdf(i,pdfJBegin,k,7) = cdiscretization_->pdf(i,pdfJBegin,k,9);
            cdiscretization_->pdf(i,pdfJBegin,k,8) = cdiscretization_->pdf(i,pdfJBegin,k,10);
            cdiscretization_->pdf(i,pdfJBegin,k,15) = cdiscretization_->pdf(i,pdfJBegin,k,17);
            cdiscretization_->pdf(i,pdfJBegin,k,16) = cdiscretization_->pdf(i,pdfJBegin,k,18);         

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
            cdiscretization_->pdf(i,pdfJEnd,k,4) = cdiscretization_->pdf(i,pdfJEnd,k,3);
            cdiscretization_->pdf(i,pdfJEnd,k,9) = cdiscretization_->pdf(i,pdfJEnd,k,7);
            cdiscretization_->pdf(i,pdfJEnd,k,10) = cdiscretization_->pdf(i,pdfJEnd,k,8); 
            cdiscretization_->pdf(i,pdfJEnd,k,17) = cdiscretization_->pdf(i,pdfJEnd,k,15);
            cdiscretization_->pdf(i,pdfJEnd,k,18) = cdiscretization_->pdf(i,pdfJEnd,k,16);
        }
    }

//set bottom left edge
for(int j = pdfJBegin; j<pdfJEnd+1; j++)
{                                                         
                cdiscretization_->pdf(pdfIBegin,j,pdfKBegin,7) = cdiscretization_->pdf(pdfIBegin,j,pdfKBegin,9);
                cdiscretization_->pdf(pdfIBegin,j,pdfKBegin,10) = cdiscretization_->pdf(pdfIBegin,j,pdfKBegin,8);
                cdiscretization_->pdf(pdfIBegin,j,pdfKBegin,15) = cdiscretization_->pdf(pdfIBegin,j,pdfKBegin,17);
                cdiscretization_->pdf(pdfIBegin,j,pdfKBegin,18) = cdiscretization_->pdf(pdfIBegin,j,pdfKBegin,16);
}
//set top left edge
for(int j = pdfJBegin+1; j<pdfJEnd; j++)
{                                                         
                cdiscretization_->pdf(pdfIBegin,j,pdfKEnd,7) = cdiscretization_->pdf(pdfIBegin,j,pdfKEnd,9);
                cdiscretization_->pdf(pdfIBegin,j,pdfKEnd,10) = cdiscretization_->pdf(pdfIBegin,j,pdfKEnd,8);
                cdiscretization_->pdf(pdfIBegin,j,pdfKEnd,17) = cdiscretization_->pdf(pdfIBegin,j,pdfKEnd,15);
                cdiscretization_->pdf(pdfIBegin,j,pdfKEnd,16) = cdiscretization_->pdf(pdfIBegin,j,pdfKEnd,18);
}
//set top right edge
for(int j = pdfJBegin+1; j<pdfJEnd; j++)
{

    cdiscretization_->pdf(pdfIEnd,j,pdfKEnd,9) = cdiscretization_->pdf(pdfIEnd,j,pdfKEnd,7);
    cdiscretization_->pdf(pdfIEnd,j,pdfKEnd,8) = cdiscretization_->pdf(pdfIEnd,j,pdfKEnd,10);
    cdiscretization_->pdf(pdfIEnd,j,pdfKEnd,17) = cdiscretization_->pdf(pdfIEnd,j,pdfKEnd,15);
    cdiscretization_->pdf(pdfIEnd,j,pdfKEnd,16) = cdiscretization_->pdf(pdfIEnd,j,pdfKEnd,18);
}
//
//set bottom right edge
for(int j = pdfJBegin+1; j<pdfJEnd; j++)
{
    cdiscretization_->pdf(pdfIEnd,j,pdfKBegin,9) = cdiscretization_->pdf(pdfIEnd,j,pdfKBegin,7);
    cdiscretization_->pdf(pdfIEnd,j,pdfKBegin,8) = cdiscretization_->pdf(pdfIEnd,j,pdfKBegin,10);
    cdiscretization_->pdf(pdfIEnd,j,pdfKBegin,15) = cdiscretization_->pdf(pdfIEnd,j,pdfKBegin,17);
    cdiscretization_->pdf(pdfIEnd,j,pdfKBegin,18) = cdiscretization_->pdf(pdfIEnd,j,pdfKBegin,16);
}
//set front bottom edge
for(int i = pdfIBegin+1; i<pdfIEnd;i++)
{
    double N_x = 1./4.*(cdiscretization_->pdf(i,pdfJBegin,pdfKBegin,1)-cdiscretization_->pdf(i,pdfJBegin,pdfKBegin,2));
    cdiscretization_->pdf(i,pdfJBegin,pdfKBegin,7) = cdiscretization_->pdf(i,pdfJBegin,pdfKBegin,9);
    cdiscretization_->pdf(i,pdfJBegin,pdfKBegin,8) =     cdiscretization_->pdf(i,pdfJBegin,pdfKBegin,10);
    cdiscretization_->pdf(i,pdfJBegin,pdfKBegin,11) =    cdiscretization_->pdf(i,pdfJBegin,pdfKBegin,13);
    cdiscretization_->pdf(i,pdfJBegin,pdfKBegin,14) =    cdiscretization_->pdf(i,pdfJBegin,pdfKBegin,12);
}
//set front top edge
for(int i = pdfIBegin+1; i<pdfIEnd;i++)
{
    double N_x = 1./4.*(cdiscretization_->pdf(i,pdfJBegin,pdfKEnd,1)-cdiscretization_->pdf(i,pdfJBegin,pdfKEnd,2));
    cdiscretization_->pdf(i,pdfJBegin,pdfKEnd,7) = cdiscretization_->pdf(i,pdfJBegin,pdfKEnd,9);
    cdiscretization_->pdf(i,pdfJBegin,pdfKEnd,8) =     cdiscretization_->pdf(i,pdfJBegin,pdfKEnd,10);
    cdiscretization_->pdf(i,pdfJBegin,pdfKEnd,13) =    cdiscretization_->pdf(i,pdfJBegin,pdfKEnd,11);
    cdiscretization_->pdf(i,pdfJBegin,pdfKEnd,12) =    cdiscretization_->pdf(i,pdfJBegin,pdfKEnd,14);
}
//set back bottom edge
for(int i = pdfIBegin+1; i<pdfIEnd;i++)
{
    double N_x = 1./4.*(cdiscretization_->pdf(i,pdfJEnd,pdfKBegin,1)-cdiscretization_->pdf(i,pdfJEnd,pdfKBegin,2));
    cdiscretization_->pdf(i,pdfJEnd,pdfKBegin,9) = cdiscretization_->pdf(i,pdfJEnd,pdfKBegin,7);
    cdiscretization_->pdf(i,pdfJEnd,pdfKBegin,10) =     cdiscretization_->pdf(i,pdfJEnd,pdfKBegin,8);
    cdiscretization_->pdf(i,pdfJEnd,pdfKBegin,11) =    cdiscretization_->pdf(i,pdfJEnd,pdfKBegin,13);
    cdiscretization_->pdf(i,pdfJEnd,pdfKBegin,14) =    cdiscretization_->pdf(i,pdfJEnd,pdfKBegin,12);
}
//set back top edge
for(int i = pdfIBegin+1; i<pdfIEnd;i++)
{
    double N_x = 1./4.*(cdiscretization_->pdf(i,pdfJEnd,pdfKEnd,1)-cdiscretization_->pdf(i,pdfJEnd,pdfKEnd,2));
    cdiscretization_->pdf(i,pdfJEnd,pdfKEnd,9) = cdiscretization_->pdf(i,pdfJEnd,pdfKEnd,7);
    cdiscretization_->pdf(i,pdfJEnd,pdfKEnd,10) =     cdiscretization_->pdf(i,pdfJEnd,pdfKEnd,8);
    cdiscretization_->pdf(i,pdfJEnd,pdfKEnd,13) =    cdiscretization_->pdf(i,pdfJEnd,pdfKEnd,11);
    cdiscretization_->pdf(i,pdfJEnd,pdfKEnd,12) =    cdiscretization_->pdf(i,pdfJEnd,pdfKEnd,14);
}
//set front left edge
for(int k = pdfKBegin; k<pdfKEnd+1;k++)
{
    double N_z = 1./4.*(cdiscretization_->pdf(pdfIBegin,pdfJBegin,k,5)-cdiscretization_->pdf(pdfIBegin,pdfJBegin,k,6));
    cdiscretization_->pdf(pdfIBegin,pdfJBegin,k,11) =  cdiscretization_->pdf(pdfIBegin,pdfJBegin,k,13);
    cdiscretization_->pdf(pdfIBegin,pdfJBegin,k,12) = cdiscretization_->pdf(pdfIBegin,pdfJBegin,k,14);
    cdiscretization_->pdf(pdfIBegin,pdfJBegin,k,15) = cdiscretization_->pdf(pdfIBegin,pdfJBegin,k,17);
    cdiscretization_->pdf(pdfIBegin,pdfJBegin,k,16) = cdiscretization_->pdf(pdfIBegin,pdfJBegin,k,18);
}
//set back left edge
for(int k = pdfKBegin; k<pdfKEnd+1;k++)
{
    double N_z = 1./4.*(cdiscretization_->pdf(pdfIBegin,pdfJEnd,k,5)-cdiscretization_->pdf(pdfIBegin,pdfJEnd,k,6));
    cdiscretization_->pdf(pdfIBegin,pdfJEnd,k,11) =  cdiscretization_->pdf(pdfIBegin,pdfJEnd,k,13);
    cdiscretization_->pdf(pdfIBegin,pdfJEnd,k,12) = cdiscretization_->pdf(pdfIBegin,pdfJEnd,k,14);
    cdiscretization_->pdf(pdfIBegin,pdfJEnd,k,17) = cdiscretization_->pdf(pdfIBegin,pdfJEnd,k,15);
    cdiscretization_->pdf(pdfIBegin,pdfJEnd,k,18) = cdiscretization_->pdf(pdfIBegin,pdfJEnd,k,16);
}

//set front right edge
for(int k = pdfKBegin; k<pdfKEnd+1;k++)
{
    double N_z = 1./4.*(cdiscretization_->pdf(pdfIEnd,pdfJBegin,k,5)-cdiscretization_->pdf(pdfIEnd,pdfJBegin,k,6));
    cdiscretization_->pdf(pdfIEnd,pdfJBegin,k,13) =  cdiscretization_->pdf(pdfIEnd,pdfJBegin,k,11);
    cdiscretization_->pdf(pdfIEnd,pdfJBegin,k,14) = cdiscretization_->pdf(pdfIEnd,pdfJBegin,k,12);
    cdiscretization_->pdf(pdfIEnd,pdfJBegin,k,15) = cdiscretization_->pdf(pdfIEnd,pdfJBegin,k,17);
    cdiscretization_->pdf(pdfIEnd,pdfJBegin,k,16) = cdiscretization_->pdf(pdfIEnd,pdfJBegin,k,18);
}
//set front right edge
for(int k = pdfKBegin; k<pdfKEnd+1;k++)
{
    double N_z = 1./4.*(cdiscretization_->pdf(pdfIEnd,pdfJEnd,k,5)-cdiscretization_->pdf(pdfIEnd,pdfJEnd,k,6));
    cdiscretization_->pdf(pdfIEnd,pdfJEnd,k,13) =  cdiscretization_->pdf(pdfIEnd,pdfJEnd,k,11);
    cdiscretization_->pdf(pdfIEnd,pdfJEnd,k,14) = cdiscretization_->pdf(pdfIEnd,pdfJEnd,k,12);
    cdiscretization_->pdf(pdfIEnd,pdfJEnd,k,17) = cdiscretization_->pdf(pdfIEnd,pdfJEnd,k,15);
    cdiscretization_->pdf(pdfIEnd,pdfJEnd,k,18) = cdiscretization_->pdf(pdfIEnd,pdfJEnd,k,16);
}
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

            for (int k = cdiscretization_->pdfKBegin(); k < cdiscretization_->pdfKEnd(); k++)
            {
                            double p = 0.;
                            double u = 0.;
                            double v = 0.;
                            double w = 0.;
                            double rho = 0.;
                for (int l = 0; l < 19; l++)
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
