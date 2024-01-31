#include "discretization/central_grid.h"
#include "central_grid.h"

CentralGrid::CentralGrid(std::array<int, 2> nCells,
                  std::array<double, 2> meshWidth):
                  nCells_(nCells),meshWidth_(meshWidth),
                    p_({nCells_[0], nCells_[1]}, {0.,0.}, meshWidth),
                    u_({nCells_[0], nCells_[1]}, {0.,0.}, meshWidth),
                    v_({nCells_[0], nCells_[1]}, {0.,0.}, meshWidth),
                    rho_({nCells_[0], nCells_[1]}, {0.,0.}, meshWidth),
                    pdf_({nCells_[0], nCells_[1], 9}, {0.,0.}, meshWidth),
                    pdfold_({nCells_[0], nCells_[1], 9}, {0.,0.}, meshWidth),
                    pdfeq_({nCells_[0], nCells_[1], 9}, {0.,0.}, meshWidth)
{
        cix_[0] = 0.;
        ciy_[0] = 0.; 
        wi_[0] = 4./9.;

        cix_[1] = 1.;
        ciy_[1] = 0.;
        
        cix_[2] = 0.;
        ciy_[2] = 1.;
        
        cix_[3] = -1.;
        ciy_[3] = 0.;
        
        cix_[4] = 0.;
        ciy_[4] = -1.;
        
        cix_[5] = 0.70710678118;
        ciy_[5] = 0.70710678118;
        
        cix_[6] = -0.70710678118;
        ciy_[6] = 0.70710678118;        
        
        
        cix_[7] = -0.70710678118;
        ciy_[7] = -0.70710678118;         
        
        cix_[8] = 0.70710678118;
        ciy_[8] = -0.70710678118;


        for(int i = 1; i < 5; i++)
        {
            wi_[i]   = 1./9.;
        }

        for(int i = 5; i<9;i++)
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