#pragma once 
#include "output_writer/output_writer.h"

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <cstdlib>
#include <iostream>
#include <memory>

class OutputWriterParaview : public OutputWriter
{ 
    public:
        OutputWriterParaview(std::shared_ptr<Discretization> discretization);
    
    
        void writeFile(double currentTime);

    private:
        vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter_;
};


