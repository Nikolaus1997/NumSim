#pragma once

#include "storage/pdffield.h"
#include "storage/fieldvariable.h"
#include <cmath>
class CentralGrid
{
    public:
        CentralGrid(std::array<int, 3> nCells,
                  std::array<double, 3> meshWidth);
        
        int pdfJEnd() const;
        int pdfJBegin() const;
        int pdfIEnd() const;
        int pdfIBegin() const;
        int pdfKEnd() const;
        int pdfKBegin() const;        
        
        int pdfeqJEnd() const;
        int pdfeqJBegin() const;
        int pdfeqIEnd() const;
        int pdfeqIBegin() const;
        int pdfeqKEnd() const;
        int pdfeqKBegin() const;

        int pJEnd() const;
        int pJBegin() const;
        int pIEnd() const;
        int pIBegin() const; 
        int pKEnd() const;
        int pKBegin() const;         
        
        int uJEnd() const;
        int uJBegin() const;
        int uIEnd() const;
        int uIBegin() const; 
        int uKEnd() const;
        int uKBegin() const;               
        
        int vJEnd() const;
        int vJBegin() const;
        int vIEnd() const;
        int vIBegin() const;  
        int vKEnd() const;
        int vKBegin() const;         

        int wJEnd() const;
        int wJBegin() const;
        int wIEnd() const;
        int wIBegin() const;   
        int wKEnd() const;
        int wKBegin() const;                               
        
        int rhoJEnd() const;
        int rhoJBegin() const;
        int rhoIEnd() const;
        int rhoIBegin() const;
        int rhoKEnd() const;
        int rhoKBegin() const;        

        double ci_x(int i) const;
        double ci_y(int i) const; 
        double ci_z(int i) const;         
        double & ci_x(int i) ;
        double & ci_y(int i) ;  
        double & ci_z(int i) ;                
        double wi(int i) const;        
        double v(int i, int j, int k) const;
        double u(int i, int j, int k) const;
        double w(int i, int j, int k) const;        
        double p(int i, int j, int k) const;
        double rho(int i, int j, int k) const;
        double pdf(int i, int j, int k, int l) const;
        double pdfold(int i, int j, int k, int l) const;
        double pdfeq(int i, int j, int k, int l) const;

        double & v(int i, int j, int k);
        double & u(int i, int j, int k);
        double & w(int i, int j, int k);        
        double & p(int i, int j, int k);
        double & rho(int i, int j, int k);
        double & pdf(int i, int j, int k, int l);        
        double & pdfold(int i, int j, int k, int l);          
        double & pdfeq(int i, int j, int k, int l);        

        double dx() const;
        double dy() const;
        const std::array<double, 3> meshWidth() const;
        const std::array<int, 3> nCells() const;

        const FieldVariable & rho() const;
        const FieldVariable & u() const;
        const FieldVariable & v() const;
        const FieldVariable & w() const;        
        const FieldVariable & p() const;
        const PdfField      & pdf() const;
        const PdfField      & pdfold() const;        
        const PdfField      & pdfeq() const;    

        const std::array<double, 3> meshWidth_;
        const std::array<int, 3> nCells_;
    
        std::array<double,19> cix_;
        std::array<double,19> ciy_;
        std::array<double,19> ciz_;
        std::array<double,19> wi_;

        FieldVariable p_;
        FieldVariable u_;
        FieldVariable v_;
        FieldVariable w_;
        FieldVariable rho_;

        PdfField    pdfeq_;
        PdfField    pdf_;
        PdfField    pdfold_;        
};