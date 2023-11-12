#pragma once
#include <memory>
#include "discretization/discretization.h"


class OutputWriter
{
    public:
        OutputWriter(std::shared_ptr<Discretization> discretization);

        virtual void writeFile(double currentTime);
    
    protected:
        std::shared_ptr< Discretization > 	discretization_;
        int 	                                    fileNo_;
};