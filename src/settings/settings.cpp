#include "settings/settings.h"
#include <fstream>   // for file operations
#include <iomanip>


void Settings::loadFromFile(std::string filename)
{
  std::string parameterName;
  std::string parameterValue;

  // open file
  std::ifstream file(filename.c_str(), std::ios::in);

  // check if file is open
  if (!file.is_open())
  {
    std::cout << "Could not open parameter file \"" << filename << "\"." << std::endl;
    return;
  }

  // loop over lines of file
  for (int lineNo = 0;; lineNo++)
  {
    // read line
    std::string line;
    getline(file, line);

    // at the end of the file break for loop
    if (file.eof())
      break;
    //skip lines that start with #
    if (line[0]=='#')
    {
      continue;
    }

    //skip lines with no = in them
    if (line.find('=')==std::string::npos)
    {
      continue;
    }
    
    //erase blank spaces at the beginning
    while (line.find_first_of("\t")==0||line.find_first_of(" ")==0)
    {
      line.erase(0,1);
    }
    
    parameterName = line.substr(0,line.find_first_of('='));
    
    if (parameterName.find_first_of(" \t") != std::string::npos)
    {
      parameterName.erase(parameterName.find_first_of(" \t"));
    }
    
    
    parameterValue = line.substr(line.find_first_of('=')+1);
    
    while(parameterValue.find(' ')!= std::string::npos)
    {
      parameterValue.erase(parameterValue.find_first_of(' '),1);
    }

    if(parameterValue.find('#')!= std::string::npos)
    {
      parameterValue.erase(parameterValue.find('#'));
    }

    if (parameterName=="endTime")
    {
      endTime = atoi(parameterValue.c_str());
    }

    if (parameterName=="physicalSizeX")
    {
      physicalSize[0] = atoi(parameterValue.c_str());
    }

    if (parameterName=="physicalSizeY")
    {
      physicalSize[1] = atoi(parameterValue.c_str());
    }
    
    if (parameterName=="physicalSizeZ")
    {
      physicalSize[2] = atoi(parameterValue.c_str());
    }
    if (parameterName=="re")
    {
      re = atoi(parameterValue.c_str());
    }    

    if (parameterName=="gX")
    {
      g[0] = atoi(parameterValue.c_str());
    }

    if (parameterName=="gY")
    {
     g[1] = atoi(parameterValue.c_str());
    }

    if (parameterName=="dirichletBottomX")
    {
     dirichletBcBottom[0] = atof(parameterValue.c_str());  
    }

    if (parameterName=="dirichletBottomY")
    {
      dirichletBcBottom[1] = atof(parameterValue.c_str());  
    }    

    if (parameterName=="dirichletBottomZ")
    {
      dirichletBcBottom[2] = atof(parameterValue.c_str());  
    }

    if (parameterName=="dirichletTopX")
    {
        dirichletBcTop[0] = atof(parameterValue.c_str());  
    }

    if (parameterName=="dirichletTopY")
    {
        dirichletBcTop[1] = atof(parameterValue.c_str());      
    }

    if (parameterName=="dirichletTopZ")
    {
        dirichletBcTop[2] = atof(parameterValue.c_str());      
    }

    if (parameterName=="dirichletLeftX")
    {
      dirichletBcLeft[0] = atof(parameterValue.c_str());
    }

    if (parameterName=="dirichletLeftY")
    {
        dirichletBcLeft[1] = atof(parameterValue.c_str());
    }

    if (parameterName=="dirichletLeftZ")
    {
        dirichletBcLeft[2] = atof(parameterValue.c_str());
    }

    if (parameterName=="dirichletRightX")
    {
       dirichletBcRight[0] = atof(parameterValue.c_str()); 
    }

    if (parameterName=="dirichletRightY")
    {
       dirichletBcRight[1] = atof(parameterValue.c_str()); 
    }

    if (parameterName=="dirichletRightZ")
    {
       dirichletBcRight[2] = atof(parameterValue.c_str()); 
    }    

    if (parameterName=="nCellsX")
    {
         nCells[0] = atof(parameterValue.c_str());
    }

    if (parameterName=="nCellsY")
    {
         nCells[1] = atof(parameterValue.c_str());
    }

    if (parameterName=="nCellsZ")
    {
         nCells[2] = atof(parameterValue.c_str());
    }

    if (parameterName=="L_lbm")
    {
         L_lbm = atof(parameterValue.c_str());
    }    

    if (parameterName=="useDonorCell")
    {
      if(parameterValue=="TRUE"||"true"){
        useDonorCell = true;
      }
      else
      {
        useDonorCell = false;
      }
    }

    if (parameterName=="tau")
    {
        tau = atof(parameterValue.c_str());
    }

    if (parameterName=="alpha")
    {
        alpha = atof(parameterValue.c_str());
    }

    if (parameterName=="maximumDt")
    {
        maximumDt = atof(parameterValue.c_str());
    }

    if (parameterName=="pressureSolver")
    {
        pressureSolver = parameterValue;
    }

    if (parameterName=="omega")
    {
        omega = atof(parameterValue.c_str());
    }

    if (parameterName=="epsilon")
    {
        epsilon = atof(parameterValue.c_str());
    }

    if (parameterName=="maximumNumberOfIterations")
    {
       maximumNumberOfIterations = atof(parameterValue.c_str());
    }

    // print line
    std::cout << "line " << lineNo << ": " << line << std::endl;
  }
}

void Settings::printSettings()
{
  std::cout << "Settings: " << std::endl
    << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl
    << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0] << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt << std::endl
    << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1]  << ")"
    << ", top: ("  << dirichletBcTop[0] << "," << dirichletBcTop[1]  << ")"
    << ", left: ("  << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
    << ", right: ("  << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl
    << "  useDonorCell: " << std::boolalpha << useDonorCell << ", alpha: " << alpha << std::endl
    << "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}