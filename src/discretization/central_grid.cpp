#include "discretization/central_grid.h"
#include "central_grid.h"

CentralGrid::CentralGrid(std::array<int, 2> nCells,
                  std::array<double, 2> meshWidth):
                  nCells_(nCells),meshWidth_(meshWidth),
                    p_({nCells_[0], nCells_[1]}, {0.,0.}, meshWidth),
                    u_({nCells_[0], nCells_[1]}, {0.,0.}, meshWidth),
                    v_({nCells_[0], nCells_[1]}, {0.,0.}, meshWidth),
                    w_({nCells_[0], nCells_[1]}, {0.,0.}, meshWidth),
                    rho_({nCells_[0], nCells_[1]}, {0.,0.}, meshWidth),
                    pdf_({nCells_[0], nCells_[1], 19}, {0.,0.}, meshWidth),
                    pdfold_({nCells_[0], nCells_[1], 19}, {0.,0.}, meshWidth),
                    pdfeq_({nCells_[0], nCells_[1], 19}, {0.,0.}, meshWidth)
{
        cix_[0] = 0.;
        ciy_[0] = 0.; 
        ciz_[0] = 0.; 
        wi_[0] = 1./3.;

        cix_[1] = 1.;
        ciy_[1] = 0.;
        ciz_[1] = 0.;

        cix_[2] = -1.;
        ciy_[2] = 0.;
        ciz_[2] = 0.;        
        
        cix_[3] = 0.;
        ciy_[3] = 1.;
        ciz_[3] = 0.;

        cix_[4] = 0.;
        ciy_[4] = -1.;
        ciz_[4] = 0.;

        cix_[5] = 0.;//0.70710678118;
        ciy_[5] = 0.;//0.70710678118;
        ciz_[5] = 1.;        

        cix_[6] = 0.;//0.70710678118;
        ciy_[6] = 0.;//0.70710678118;        
        ciz_[6] = -1.;        
        
        cix_[7] = 1.;//0.70710678118;
        ciy_[7] = 1.;//0.70710678118;         
        ciz_[7] = 0.;        

        cix_[8] = -1.;//0.70710678118;
        ciy_[8] = 1.;//0.70710678118;
        ciz_[8] = 0.;

        cix_[9] = -1.;
        ciy_[9] = -1.;
        ciz_[9] = 0.;

        cix_[10] = 1.;
        ciy_[10] = -1.;
        ciz_[10] = 0.;

        cix_[11] = 1.;
        ciy_[11] = 0.;
        ciz_[11] = 1.;

        cix_[12] = 1.;
        ciy_[12] = 0.;
        ciz_[12] = -1.;

        cix_[13] = -1.;
        ciy_[13] = 0.;
        ciz_[13] = -1.;        

        cix_[14] = -1.;
        ciy_[14] = 0.;
        ciz_[14] = 1.;

        cix_[15] = 0.;
        ciy_[15] = 1.;
        ciz_[15] = 1.; 

        cix_[16] = 0.;
        ciy_[16] = 1.;
        ciz_[16] = -1.; 

        cix_[17] = 0.;
        ciy_[17] = -1.;
        ciz_[17] = -1.;                

        cix_[18] = 0.;
        ciy_[18] = -1.;
        ciz_[18] = 1.; 

        for(int i = 1; i < 7; i++)
        {
            wi_[i]   = 1./18.;
        }

        for(int i = 8; i<19;i++)
        {
            wi_[i] = 1./36.;
        }
};

/**
 * get mesh width in x-direction
 * @return mesh width in x-direction
 */
double CentralGrid::dx() const
{
    return meshWidth_[0];
};

/**
 * get mesh width in y-direction
 * @return mesh width in y-direction
 */
double CentralGrid::dy() const
{
    return meshWidth_[1];
};


const std::array<double, 2> CentralGrid::meshWidth() const
{
    return meshWidth_;
};

const std::array<int, 2> CentralGrid::nCells() const
{
    return nCells_;
};

int CentralGrid::pdfIBegin() const
{
    return 0;
};

int CentralGrid::pdfIEnd() const
{
    return nCells_[0];
};

int CentralGrid::pdfJBegin() const
{
    return 0;
};

int CentralGrid::pdfJEnd() const
{
    return nCells_[1];
};

const PdfField &CentralGrid::pdf() const
{
    return pdf_;
};
const PdfField &CentralGrid::pdfold() const
{
    return pdfold_;
};
int CentralGrid::pdfeqIBegin() const
{
    return 0;
};

int CentralGrid::pdfeqIEnd() const
{
    return nCells_[0];
};

int CentralGrid::pdfeqJBegin() const
{
    return 0;
};

int CentralGrid::pdfeqJEnd() const
{
    return nCells_[1];
};

const PdfField &CentralGrid::pdfeq() const
{
    return pdfeq_;
};

int CentralGrid::pIBegin() const
{
    return 0;
};

int CentralGrid::pIEnd() const
{
    return nCells_[0];
};

int CentralGrid::pJBegin() const
{
    return 0;
};

int CentralGrid::pJEnd() const
{
    return nCells_[1];
};


const FieldVariable &CentralGrid::p() const
{
    return p_;
};

int CentralGrid::rhoIBegin() const
{
    return 0;
};

int CentralGrid::rhoIEnd() const
{
    return nCells_[0];
};

int CentralGrid::rhoJBegin() const
{
    return 0;
};

int CentralGrid::rhoJEnd() const
{
    return nCells_[1];
};


const FieldVariable &CentralGrid::rho() const
{
    return rho_;
};

int CentralGrid::uIBegin() const
{
    return 0;
};

int CentralGrid::uIEnd() const
{
    return nCells_[0];
};

int CentralGrid::uJBegin() const
{
    return 0;
};

int CentralGrid::uJEnd() const
{
    return nCells_[1];
};


const FieldVariable &CentralGrid::u() const
{
    return u_;
};

int CentralGrid::vIBegin() const
{
    return 0;
};

int CentralGrid::vIEnd() const
{
    return nCells_[0];
};

double CentralGrid::v(int i, int j) const
{
    return v_(i-vIBegin(),j-vJBegin());
}

double CentralGrid::ci_x(int i) const
{
    return cix_[i];
}
double CentralGrid::ci_y(int i) const
{
    return ciy_[i];
}
double &CentralGrid::ci_x(int i) 
{
    return cix_[i];
}
double &CentralGrid::ci_y(int i) 
{
    return ciy_[i];
}
double CentralGrid::wi(int i) const
{
    return wi_[i];
}

double CentralGrid::p(int i, int j) const
{
    return p_(i-pIBegin(),j-pJBegin());
}

double CentralGrid::u(int i, int j) const
{
    return u_(i-uIBegin(),j-uJBegin());
}

double CentralGrid::rho(int i, int j) const
{
    return rho_(i-rhoIBegin(),j-rhoJBegin());
}

double CentralGrid::pdf(int i, int j, int k) const
{
    return pdf_(i-pdfIBegin(),j-pdfJBegin(), k);
}

double CentralGrid::pdfold(int i, int j, int k) const
{
    return pdfold_(i-pdfIBegin(),j-pdfJBegin(), k);
}

double CentralGrid::pdfeq(int i, int j, int k) const
{
    return pdfeq_(i-pdfeqIBegin(),j-pdfeqJBegin(), k);
}

int CentralGrid::vJBegin() const
{
    return 0;
};

int CentralGrid::vJEnd() const
{
    return nCells_[1];
};

const FieldVariable &CentralGrid::v() const
{
    return v_;
};


double &CentralGrid::p(int i, int j)
{
    return p_(i - pIBegin(), j - pJBegin());
};

double &CentralGrid::u(int i, int j)
{
    return u_(i - uIBegin(), j - uJBegin());
};

double &CentralGrid::v(int i, int j)
{
    return v_(i - vIBegin(), j - vJBegin());
};

double &CentralGrid::rho(int i, int j)
{
    return rho_(i - rhoIBegin(), j - rhoJBegin());
};

double &CentralGrid::pdf(int i, int j, int k)
{
    return pdf_(i - pdfIBegin(), j - pdfJBegin(), k);
};
double &CentralGrid::pdfold(int i, int j, int k)
{
    return pdfold_(i - pdfIBegin(), j - pdfJBegin(), k);
};

double &CentralGrid::pdfeq(int i, int j, int k)
{
    return pdfeq_(i - pdfeqIBegin(), j - pdfeqJBegin(), k);
};