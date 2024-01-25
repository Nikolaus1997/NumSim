#pragma once

#include "storage/pdffield.h"
#include "storage/fieldvariable.h"

class CentralGrid
{
    public:
        CentralGrid(std::array<int, 2> nCells,
                  std::array<double, 2> meshWidth);
        int pdfJEnd() const;
        int pdfJBegin() const;
        int pdfIEnd() const;
        int pdfIBegin() const;
        const FieldVariable & u() const;
        const FieldVariable & v() const;
        const FieldVariable & p() const;
        const PdfField      & pdf() const;
    protected:
        const std::array<double, 2> meshWidth_;
        const std::array<int, 2> nCells_;
        FieldVariable p_;
        FieldVariable u_;
        FieldVariable v_;
        PdfField    pdf_;
};