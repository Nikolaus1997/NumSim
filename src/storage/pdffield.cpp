#include "storage/pdffield.h"

PdfField::PdfField(std::array<int, 3> size,
                             std::array<double, 2> origin,
                             std::array<double, 2> meshWidth) :
    Array3D(size),
    origin_(origin),
    meshWidth_(meshWidth)
{
};