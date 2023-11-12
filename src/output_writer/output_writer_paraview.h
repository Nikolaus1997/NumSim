#pragma once 
#include "output_writer/output_writer.h"

#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

class OutputWriterParaview : public OutputWriter
{ 
    public:
        OutputWriterParaview(std::shared_ptr<Discretization> discretization);
    
    
        void writeFile(double currentTime);

    private:
        vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter_;
};


