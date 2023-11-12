#pragma once
#include "output_writer/output_writer.h"

class OutputWriterText: public OutputWriter
{
    public:
        void writeFile(double currentTime);

        void writePressureFile();
};