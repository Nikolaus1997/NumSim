#include "output_writer_paraview.h"

#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

OutputWriterParaview::OutputWriterParaview(std::shared_ptr<Discretization> discretization)
:OutputWriter(discretization)
{
    // Create a vtkWriter_
    vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
}
//TODO