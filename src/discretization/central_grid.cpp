#include "discretization/central_grid.h"

CentralGrid::CentralGrid(std::array<int, 2> nCells,
                  std::array<double, 2> meshWidth, std::array<double,2> NumberLbmDirections):
                  nCells_(nCells),meshWidth_(meshWidth),
                    p_({nCells_[0], nCells_[1]}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
                    u_({nCells_[0], nCells_[1]}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
                    v_({nCells_[0], nCells_[1]}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
                    pdf_({nCells_[0]+2, nCells_[1]+2, 9}, {meshWidth[0] / 2., meshWidth[1] / 2.}, meshWidth),
                    ci_(NumberLbmDirections)
{
    for(int i = 0; i < 9; i++)
    {

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

const FieldVariable &CentralGrid::p() const
{
    return p_;
};

const FieldVariable &CentralGrid::u() const
{
    return u_;
};

const FieldVariable &CentralGrid::v() const
{
    return v_;
};