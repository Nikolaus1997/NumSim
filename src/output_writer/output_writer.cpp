#include "output_writer.h"

OutputWriter::OutputWriter(std::shared_ptr<Discretization> discretization)
:discretization_(discretization), fileNo_(-1)
{
    // create "out" subdirectory if it does not yet exist
    int returnValue = system("mkdir -p out");
    if (returnValue != 0)
    std::cout << "Could not create subdirectory \"out\"." << std::endl;
}

