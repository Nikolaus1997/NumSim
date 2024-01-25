#include "discretization/central_grid.h"
#include "central_grid.h"

CentralGrid::CentralGrid(std::array<int, 2> nCells,
                  std::array<double, 2> meshWidth):
                  nCells_(nCells),meshWidth_(meshWidth),
                    p_({nCells_[0], nCells_[1]}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
                    u_({nCells_[0], nCells_[1]}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
                    v_({nCells_[0], nCells_[1]}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
                    rho_({nCells_[0], nCells_[1]}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
                    pdf_({nCells_[0]+2, nCells_[1]+2, 9}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
                    ci_({9,2}),
                    wi_({9,1})
{
    ci_(0,1) = 0.;
    ci_(0,0) = 0.; 
    wi_(0,0) = 4/9;
    for(int i = 1; i < 5; i++)
    {
        ci_(i,0) = cos((i-1)/2*M_PI);
        ci_(i,1) = sin((i-1)/2*M_PI);
        wi_(i,0)   = 1/9;
    }

    for(int i = 5; i<9;i++)
    {
        ci_(i,0) = cos((2*i-9)/2*M_PI);
        ci_(i,1) = sin((2*i-9)/2*M_PI);
        wi_(i,0) = 1/36;
    }
};



int CentralGrid::pdfIBegin() const
{
    return -1;
};

int CentralGrid::pdfIEnd() const
{
    return nCells_[0] + 1;
};

int CentralGrid::pdfJBegin() const
{
    return -1;
};

int CentralGrid::pdfJEnd() const
{
    return nCells_[1] + 1;
};

const PdfField &CentralGrid::pdf() const
{
    return pdf_;
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