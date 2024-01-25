#pragma once

#include "storage/pdffield.h"
#include "storage/fieldvariable.h"
#include <cmath>
class CentralGrid
{
    public:
        CentralGrid(std::array<int, 2> nCells,
                  std::array<double, 2> meshWidth);
        
        int pdfJEnd() const;
        int pdfJBegin() const;
        int pdfIEnd() const;
        int pdfIBegin() const;
        
        int pJEnd() const;
        int pJBegin() const;
        int pIEnd() const;
        int pIBegin() const; 
        
        int uJEnd() const;
        int uJBegin() const;
        int uIEnd() const;
        int uIBegin() const;       
        
        int vJEnd() const;
        int vJBegin() const;
        int vIEnd() const;
        int vIBegin() const;                  
        
        int rhoJEnd() const;
        int rhoJBegin() const;
        int rhoIEnd() const;
        int rhoIBegin() const;

        double v(int i, int j) const;
        double u(int i, int j) const;
        double p(int i, int j) const;
        double rho(int i, int j) const;
        double pdf(int i, int j, int k) const;

        double & v(int i, int j);
        double & u(int i, int j);
        double & p(int i, int j);
        double & rho(int i, int j);
        double & pdf(int i, int j, int k);        

        const FieldVariable & rho() const;
        const FieldVariable & u() const;
        const FieldVariable & v() const;
        const FieldVariable & p() const;
        const PdfField      & pdf() const;
    
    protected:
        const std::array<double, 2> meshWidth_;
        const std::array<int, 2> nCells_;
        
        Array2D wi_;
        Array2D ci_;
        
        FieldVariable p_;
        FieldVariable u_;
        FieldVariable v_;
        FieldVariable rho_;

        PdfField    pdf_;
};