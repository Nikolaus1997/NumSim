#include "output_writer_paraview.h"

#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

OutputWriterParaview::OutputWriterParaview(std::shared_ptr<Discretization> discretization)
:OutputWriter(discretization)
{
    // Create a vtkWriter_
    vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
};

void OutputWriterParaview::writeFile(double currentTime)
{
    // Assemble the filename
    std::stringstream fileName;
    if (vtkWriter_ != nullptr) {
    vtkWriter_->GetDefaultFileExtension();
    } else {
    std::cout << "vtkWriter_ is null" << std::endl;
    }
    fileName << "out/output_" << std::setw(4) << setfill('0')<< fileNo_  << "." << vtkWriter_->GetDefaultFileExtension();

    // increment fileNo
    fileNo_++;
    std::cout << "Write file \"" << fileName.str() << "\"." << std::endl;

    // assign the new file name to the output vtkWriter
    vtkWriter_->SetFileName(fileName.str().c_str());

    // initialize data set that will be output to the file
    vtkSmartPointer<vtkImageData> dataSet = vtkSmartPointer<vtkImageData>::New();
    dataSet->SetOrigin(0, 0, 0);

    // set spacing of mesh
    const double dx = discretization_->meshWidth()[0];
    const double dy = discretization_->meshWidth()[1];
    const double dz = 1;
    dataSet->SetSpacing(dx, dy, dz);

    // set number of points in each dimension, 1 cell in z direction
    std::array<int,2> nCells = discretization_->nCells();
    dataSet->SetDimensions(nCells[0]+1, nCells[1]+1, 1);

    // add pressure field variable
    // ---------------------------
    vtkSmartPointer<vtkDoubleArray> arrayPressure = vtkDoubleArray::New();

    // the pressure is a scalar which means the number of components is 1
    arrayPressure->SetNumberOfComponents(1);

    // Set the number of pressure values and allocate memory for it. We already know the number, it has to be the same as there are nodes in the mesh.
    arrayPressure->SetNumberOfTuples(dataSet->GetNumberOfPoints());

    arrayPressure->SetName("pressure");

    // loop over the nodes of the mesh and assign an artifical value that changes with fileNo
    int index = 0;   // index for the vtk data structure, will be incremented in the inner loop
    for (int j = 0; j < nCells[1]+1; j++)
    {
        for (int i = 0; i < nCells[0]+1; i++, index++)
        {
        const double x = i*dx;
        const double y = j*dy;

        arrayPressure->SetValue(index, discretization_->p().interpolateAt(x,y));
        }
    }

    // now, we should have added as many values as there are points in the vtk data structure
    assert(index == dataSet->GetNumberOfPoints());  
    
    // Add the field variable to the data set
    dataSet->GetPointData()->AddArray(arrayPressure);
    
    // add current time 
    vtkSmartPointer<vtkDoubleArray> arrayTime = vtkDoubleArray::New();
    arrayTime->SetName("TIME");
    arrayTime->SetNumberOfTuples(1);
    arrayTime->SetTuple1(0, currentTime);
    dataSet->GetFieldData()->AddArray(arrayTime);

    // Remove unused memory
    dataSet->Squeeze();

    // Write the data
    vtkWriter_->SetInputData(dataSet);

    //vtkWriter->SetDataModeToAscii();     // comment this in to get ascii text files: those can be checked in an editor
    vtkWriter_->SetDataModeToBinary();      // set file mode to binary files: smaller file sizes

    // finally write out the data
    vtkWriter_->Write();
}



//TODO