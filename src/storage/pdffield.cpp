#include "storage/pdffield.h"

PdfField::PdfField(std::array<int, 4> size,
                             std::array<double, 3> origin,
                             std::array<double, 3> meshWidth) :
    Array4D(size),
    origin_(origin),
    meshWidth_(meshWidth)
{
};