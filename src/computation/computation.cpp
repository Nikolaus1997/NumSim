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
    //settings_.loadFromFile(filename);
    // Print settings

    //settings_.printSettings();
 

    // Initialize discretization
    for (int i = 0; i < 2; i++)
        meshWidth_[i] = settings_.physicalSize[i] / settings_.nCells[i];

        // doLBM = true;
        cdiscretization_ = std::make_shared<LbmDiscretization>(settings_.nCells, meshWidth_);

        cs_ =(double) 1./sqrt(3);
        nu_ = settings_.dirichletBcTop[0]*settings_.L_lbm/settings_.re;
        tau_ = (6*nu_+1)/2.;
        //tau_ =0.65;
        dt_ = 1.;


//(settings_.nCells[1]/settings_.re)*settings_.dirichletBcTop[0]+.5;


        //tau_m = 1./4.*1/(tau_-0.5)+0.5;
        // tau_om_m = 1./tau_m;
        // tau_om_p = 1./tau_;
        // tau_om   = 1./tau_;
        // dt_ = meshWidth_[0]*settings_.dirichletBcTop[0];
        // tau_ = 1.25;

                std::cout<<tau_<<std::endl;
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
void Computation::runSimulation() {
    int t_iter = 0;
    double time = 0.0;
    double holder = settings_.dirichletBcTop[0];
    double dtp_ = settings_.re*nu_/(settings_.nCells[0]*settings_.nCells[0]);
    int dtwrite_ = settings_.deltawrite;
            std::cout<<dtp_<<std::endl;
    while (time < settings_.endTime){







        // settings_.dirichletBcTop[0] = 1000.*(1.0-exp(-t_iter*t_iter/(2.0*10*settings_.nCells[0]*10*settings_.nCells[0])));
        //         t_iter++;

        //std::cout<<"+++++++++++++++++++"<<std::endl;     
        if(time == 0.0 ){
            // for(int i=cdiscretization_->pIBegin(); i < cdiscretization_->pIEnd();i++)
            //     {
            //         cdiscretization_->u(i,cdiscretization_->uJBegin())=settings_.dirichletBcBottom[0];
            //         cdiscretization_->u(i,cdiscretization_->uJEnd()-1)  =settings_.dirichletBcTop[0];
            //         cdiscretization_->v(i,cdiscretization_->vJBegin())=settings_.dirichletBcBottom[1];
            //         cdiscretization_->v(i,cdiscretization_->vJEnd()-1)  =settings_.dirichletBcTop[1];
            //     } 
            // for(int j=cdiscretization_->pJBegin(); j < cdiscretization_->pJEnd();j++)
            //     {
            //         cdiscretization_->u(cdiscretization_->uJBegin(),j)=settings_.dirichletBcLeft[0];
            //         cdiscretization_->u(cdiscretization_->uJEnd()-1,j)  =settings_.dirichletBcRight[0];
            //         cdiscretization_->v(cdiscretization_->vJBegin(),j)=settings_.dirichletBcLeft[1];
            //         cdiscretization_->v(cdiscretization_->vJEnd()-1,j)  =settings_.dirichletBcRight[1];
            //     } 
                            for(int i=cdiscretization_->pIBegin(); i < cdiscretization_->pIEnd();i++)
            {      
                    for(int j=cdiscretization_->pJBegin(); j < cdiscretization_->pJEnd();j++)
                    {
                        cdiscretization_->rho(i,j) = 1.;
                    }}
                    
        } 
     
        
        Collision();
        Streaming();
        //FillVelocitiesAndPressure();          
        //FillVelocitiesAndPressure();

        // settings_.dirichletBcLeft[0] = 0.01*(1.0-exp(-t_iter*t_iter/(10*settings_.nCells[1]))); 
        t_iter++;       
        LBMapplyBoundaryValues(); 

        //LBMBounceBack();        
        FillVelocitiesAndPressure();         
        //std::cout<<"*******************"<<std::endl;


        //std::cout<<"-------------------"<<std::endl;      
        //FillVelocitiesAndPressure();        
        //std::cout<<"+++++++++++++++++++"<<std::endl;
        if((t_iter % dtwrite_) == 0 ||time == 0.0){
            //outputWriterText_->writeFile(time);        
            outputWriterParaview_->writeFile(time);

        }



        // std::cout<<"                 "<<std::endl;
        // for(int k = 0; k<9;k++){
        //     std::cout<<"  dt:  "<<dt_<<"  "<<"########pdf0/2: "<<cdiscretization_->pdf(0,4,k)<<" "<<"           ########pdf0/1: "<<cdiscretization_->pdf(1,2,k)<<" "<<"            ########pdf0/0: "<<cdiscretization_->pdf(1,3,k)<<std::endl;
        // }

        time += dtp_;


        //}
    }

};

void Computation::LBMapplyBoundaryValues() {
    //top condition
    int pdfJEnd = cdiscretization_->pdfJEnd()-1;
    int pdfJBegin = cdiscretization_->pdfJBegin();
    int pdfIBegin = cdiscretization_->pdfIBegin();
    int pdfIEnd = cdiscretization_->pdfIEnd()-1;
//left
    for(int j = cdiscretization_->pdfJBegin()+1; j < cdiscretization_->pdfJEnd()-1; j++) {

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
    for(int j = cdiscretization_->pdfJBegin()+1; j < cdiscretization_->pdfJEnd()-1; j++) {

        double rho_N = 1.0/(1.0+settings_.dirichletBcRight[0])*(cdiscretization_->pdf(pdfIEnd,j,0)+cdiscretization_->pdf(pdfIEnd,j,2)+
            cdiscretization_->pdf(pdfIEnd,j,4)+2.*(cdiscretization_->pdf(pdfIEnd,j,5)+cdiscretization_->pdf(pdfIEnd,j,1)+cdiscretization_->pdf(pdfIEnd,j,8)));
        cdiscretization_->rho(pdfIEnd,j) = rho_N; 
        cdiscretization_->pdf(pdfIEnd,j,3) = cdiscretization_->pdf(pdfIEnd,j,1)-2./3.*rho_N*settings_.dirichletBcRight[0];
        cdiscretization_->pdf(pdfIEnd,j,7) = (cdiscretization_->pdf(pdfIEnd,j,5)+1./2.*(cdiscretization_->pdf(pdfIEnd,j,2)-cdiscretization_->pdf(pdfIEnd,j,4)) 
                                                -1./6.*rho_N*settings_.dirichletBcRight[0]-1./2.*rho_N*settings_.dirichletBcRight[1]);
        cdiscretization_->pdf(pdfIEnd,j,6) = (cdiscretization_->pdf(pdfIEnd,j,8)-1./2.*(cdiscretization_->pdf(pdfIEnd,j,2)-cdiscretization_->pdf(pdfIEnd,j,4))
                                                -1./6.*rho_N*settings_.dirichletBcRight[0]+1./2.*rho_N*settings_.dirichletBcRight[1]);

        // double rho_N = 1.0;//(1.0+settings_.dirichletBcRight[0])*(cdiscretization_->pdf(pdfIEnd,j,0)+cdiscretization_->pdf(pdfIEnd,j,2)+
        //     // cdiscretization_->pdf(pdfIEnd,j,4)+2.*(cdiscretization_->pdf(pdfIEnd,j,5)+cdiscretization_->pdf(pdfIEnd,j,1)+cdiscretization_->pdf(pdfIEnd,j,8)));
        // cdiscretization_->rho(pdfIEnd,j) = 1.0; 
        // cdiscretization_->u(pdfIEnd,j) = 0.0;
        // cdiscretization_->v(pdfIEnd,j) = 0.0;
        // cdiscretization_->pdf(pdfIEnd,j,3) = cdiscretization_->pdf(pdfIEnd,j,1);//-2./3.*rho_N*settings_.dirichletBcRight[0];
        // cdiscretization_->pdf(pdfIEnd,j,7) = (cdiscretization_->pdf(pdfIEnd,j,5)+1./2.*(cdiscretization_->pdf(pdfIEnd,j,2)-cdiscretization_->pdf(pdfIEnd,j,4))); 
        //                                         //-1./6.*rho_N*settings_.dirichletBcRight[0]-1./2.*rho_N*settings_.dirichletBcRight[1]);
        // cdiscretization_->pdf(pdfIEnd,j,6) = (cdiscretization_->pdf(pdfIEnd,j,8)-1./2.*(cdiscretization_->pdf(pdfIEnd,j,2)-cdiscretization_->pdf(pdfIEnd,j,4)));
        //                                         //-1./6.*rho_N*settings_.dirichletBcRight[0]+1./2.*rho_N*settings_.dirichletBcRight[1]);
   }    
//top
    for(int i = cdiscretization_->pdfIBegin()+1; i < cdiscretization_->pdfIEnd()-1; i++) {

        double rho_N = 1.0/(1.0+settings_.dirichletBcTop[1])*(cdiscretization_->pdf(i,pdfJEnd,0)+cdiscretization_->pdf(i,pdfJEnd,1)+
            cdiscretization_->pdf(i,pdfJEnd,3)+2.0*(cdiscretization_->pdf(i,pdfJEnd,2)+cdiscretization_->pdf(i,pdfJEnd,6)+cdiscretization_->pdf(i,pdfJEnd,5)));
        cdiscretization_->rho(i,pdfJEnd) = rho_N; 
        
        cdiscretization_->pdf(i,pdfJEnd,4) = cdiscretization_->pdf(i,pdfJEnd,2)-2./3.*rho_N*settings_.dirichletBcTop[1];

        cdiscretization_->pdf(i,pdfJEnd,7) = cdiscretization_->pdf(i,pdfJEnd,5)+1./2.*(cdiscretization_->pdf(i,pdfJEnd,1)-cdiscretization_->pdf(i,pdfJEnd,3)) 
                                                -1.0/6.0*rho_N*settings_.dirichletBcTop[1]-1.0/2.0*rho_N*settings_.dirichletBcTop[0];

        cdiscretization_->pdf(i,pdfJEnd,8) = cdiscretization_->pdf(i,pdfJEnd,6)-1./2.*(cdiscretization_->pdf(i,pdfJEnd,1)-cdiscretization_->pdf(i,pdfJEnd,3))
                                                -1.0/6.0*rho_N*settings_.dirichletBcTop[1]+1./2.*rho_N*settings_.dirichletBcTop[0];

     }
//Bottom
    for(int i = cdiscretization_->pdfIBegin()+1; i < cdiscretization_->pdfIEnd()-1; i++) {

        double rho_N = 1./(1.-settings_.dirichletBcBottom[1])*(cdiscretization_->pdf(i,pdfJBegin,0)+cdiscretization_->pdf(i,pdfJBegin,1)+
            cdiscretization_->pdf(i,pdfJBegin,3)+2.*(cdiscretization_->pdf(i,pdfJBegin,4)+cdiscretization_->pdf(i,pdfJBegin,7)+cdiscretization_->pdf(i,pdfJBegin,8)));
        cdiscretization_->rho(i,pdfJBegin) = rho_N; 
        cdiscretization_->pdf(i,pdfJBegin,2) = cdiscretization_->pdf(i,pdfJBegin,4)+2./3.*rho_N*settings_.dirichletBcBottom[1];
        cdiscretization_->pdf(i,pdfJBegin,5) = (cdiscretization_->pdf(i,pdfJBegin,7)-1./2.*(cdiscretization_->pdf(i,pdfJBegin,1)-cdiscretization_->pdf(i,pdfJBegin,3)) 
                                                +1./6.*rho_N*settings_.dirichletBcBottom[1]+1./2.*rho_N*settings_.dirichletBcBottom[0]);
        cdiscretization_->pdf(i,pdfJBegin,6) = (cdiscretization_->pdf(i,pdfJBegin,8)+1./2.*(cdiscretization_->pdf(i,pdfJBegin,1)-cdiscretization_->pdf(i,pdfJBegin,3))
                                                +1./6.*rho_N*settings_.dirichletBcBottom[1]-1./2.*rho_N*settings_.dirichletBcBottom[0]);

    }


//bottom left
cdiscretization_->u(0,0) = cdiscretization_->u(1,0);
cdiscretization_->v(0,0) = cdiscretization_->v(1,0);
cdiscretization_->rho(0,0) = cdiscretization_->rho(1,0);
cdiscretization_->pdf(0,0,1) = cdiscretization_->pdf(0,0,3)+(2.0/3.0)*cdiscretization_->rho(0,0)*cdiscretization_->u(0,0);
cdiscretization_->pdf(0,0,2) = cdiscretization_->pdf(0,0,4)+(2.0/3.0)*cdiscretization_->rho(0,0)*cdiscretization_->v(0,0);
cdiscretization_->pdf(0,0,5) = cdiscretization_->pdf(0,0,7)+(1.0/6.0)*cdiscretization_->rho(0,0)*cdiscretization_->u(0,0)+(1./6.)*cdiscretization_->rho(0,0)*cdiscretization_->v(0,0);
cdiscretization_->pdf(0,0,6) = 0.0;
cdiscretization_->pdf(0,0,8) = 0.0;
cdiscretization_->pdf(0,0,0) = cdiscretization_->rho(0,0)-cdiscretization_->pdf(0,0,1)-cdiscretization_->pdf(0,0,2)-cdiscretization_->pdf(0,0,3)-cdiscretization_->pdf(0,0,4)
                                -cdiscretization_->pdf(0,0,5)-cdiscretization_->pdf(0,0,6)-cdiscretization_->pdf(0,0,7)-cdiscretization_->pdf(0,0,8);
//top left
cdiscretization_->u(0,pdfJEnd) = cdiscretization_->u(1,pdfJEnd);
cdiscretization_->v(0,pdfJEnd) = cdiscretization_->v(1,pdfJEnd);
cdiscretization_->rho(0,cdiscretization_->pdfJEnd()-1) = cdiscretization_->rho(1,cdiscretization_->pdfJEnd()-1);

cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,1) = cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,3)
                                                            +(2.0/3.0)*cdiscretization_->rho(0,cdiscretization_->pdfJEnd()-1)*cdiscretization_->u(0,pdfJEnd);

cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,4) = cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,2);
                                                            -(2.0/3.0)*cdiscretization_->rho(0,cdiscretization_->pdfJEnd()-1)*cdiscretization_->v(0,pdfJEnd);

cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,8) = cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,6)
                                                            +(1.0/6.0)*cdiscretization_->rho(0,cdiscretization_->pdfJEnd()-1)*cdiscretization_->u(0,pdfJEnd)
                                                            -(1./6.)*cdiscretization_->rho(0,cdiscretization_->pdfJEnd()-1)*cdiscretization_->v(0,pdfJEnd);

cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,5) = 0.0;
cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,7) = 0.0;

cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,0) = cdiscretization_->rho(0,cdiscretization_->pdfJEnd()-1)-cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,1)
                                -cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,2)-cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,3)-cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,4)
                                -cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,5)-cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,6)
                                -cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,7)-cdiscretization_->pdf(0,cdiscretization_->pdfJEnd()-1,8);
//top right
cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1) = cdiscretization_->rho(cdiscretization_->pdfIEnd()-2,cdiscretization_->pdfJEnd()-1);
cdiscretization_->u(pdfIEnd,pdfJEnd) = cdiscretization_->u(pdfIEnd-1,pdfJEnd);
cdiscretization_->v(pdfIEnd,pdfJEnd) = cdiscretization_->v(pdfIEnd-1,pdfJEnd);

cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,3) = cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,1)
                                                            -(2.0/3.0)*cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1)*cdiscretization_->u(pdfIEnd,pdfJEnd);

cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,4) = cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,2)
                                                            -(2.0/3.0)*cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1)*cdiscretization_->v(pdfIEnd,pdfJEnd);

cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,7) = cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1,5)
                                                            -(1.0/6.0)*cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1)*cdiscretization_->u(pdfIEnd,pdfJEnd)
                                                           -(1.0/6.0)*cdiscretization_->rho(cdiscretization_->pdfIEnd()-1,cdiscretization_->pdfJEnd()-1)*cdiscretization_->v(pdfIEnd,pdfJEnd);
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
cdiscretization_->u(pdfIEnd,pdfJBegin) = cdiscretization_->u(pdfIEnd-1,pdfJBegin);
cdiscretization_->v(pdfIEnd,pdfJBegin) = cdiscretization_->v(pdfIEnd-1,pdfJBegin);
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
                                                            -cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,2)-cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,3)
                                                            -cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,4)
                                                            -cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,5)-cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,6)
                                                            -cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,7)-cdiscretization_->pdf(cdiscretization_->pdfIEnd()-1,0,8);
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
    for(int j = cdiscretization_->pdfJBegin(); j < cdiscretization_->pdfJEnd(); j++)
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
              if(k>=4){
                u += cdiscretization_->pdf(i,j,k)*cdiscretization_->ci_x(k);//*0.70710678118;
                v += cdiscretization_->pdf(i,j,k)*cdiscretization_->ci_y(k);//*0.70710678118;
              }else{
                u += cdiscretization_->pdf(i,j,k)*cdiscretization_->ci_x(k);
                v += cdiscretization_->pdf(i,j,k)*cdiscretization_->ci_y(k);
              }              
  
            //std::cout<<"###### "<<cdiscretization_->pdf(i,j,k)<< "####rho "<<cdiscretization_->rho(i,j)<<std::endl;               
            }

             cdiscretization_->rho(i,j) = (double)rho;
             cdiscretization_->p(i,j) = (double)p;
             cdiscretization_->u(i,j) =(double) u/rho;
             cdiscretization_->v(i,j) =(double) v/rho;
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
                     cdiscretization_->pdfeq(i,j,k) = cdiscretization_->wi(k)*cdiscretization_->rho(i,j)*(1.
                                                        +(cdiscretization_->u(i,j)*cdiscretization_->ci_x(k)
                                                                                    +cdiscretization_->v(i,j)*cdiscretization_->ci_y(k))/(cs_*cs_)
                                                        +1./2.*((cdiscretization_->u(i,j)*cdiscretization_->ci_x(k)
                                                                +cdiscretization_->v(i,j)*cdiscretization_->ci_y(k))/(cs_*cs_))
                                                        *((cdiscretization_->u(i,j)*cdiscretization_->ci_x(k)
                                                                +cdiscretization_->v(i,j)*cdiscretization_->ci_y(k))/(cs_*cs_))
                                                        -1./2.*(cdiscretization_->u(i,j)*cdiscretization_->u(i,j)
                                                                +cdiscretization_->v(i,j)*cdiscretization_->v(i,j))/(cs_*cs_));
                //std::cout<<"###### "<<cdiscretization_->pdfeq(i,j,k)<< "####rho "<<cdiscretization_->rho(i,j)<<std::endl;

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
                cdiscretization_->pdf(i,j,k) = cdiscretization_->pdf(i,j,k)-dt_/tau_*(cdiscretization_->pdf(i,j,k)-cdiscretization_->pdfeq(i,j,k));                
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

