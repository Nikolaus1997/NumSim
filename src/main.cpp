#include "output_writer/write_paraview_output.h"
#include "discretization/central_differences.h"
#include "settings/settings.h" 
#include "pressure_solver/pressure_solver.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/sor.h"
#include "computation/computation.h"
#include <iostream>
#include <cstdlib>

int main(int argc, char *argv[])
{
  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }
 std::cout <<R"(=========================================================================================================================)"<<std::endl;
 std::cout <<R"( ____  __. ___.____    .____     _____________________     .____    ________  ___________ ______________________________ )"<<std::endl;
 std::cout <<R"(|    |/  /|   |    |   |    |    \_   _____/\______   \    |    |   \_____  \ \_   _____//   _____/\_   _____/\______   \)"<<std::endl;
 std::cout <<R"(|      <  |   |    |   |    |     |    __)_  |       _/    |    |    /   |   \ |    __)_ \_____  \  |    __)_  |       _/)"<<std::endl;
 std::cout <<R"(|    |  \ |   |    |___|    |___  |        \ |    |   \    |    |___/    |    \|        \/        \ |        \ |    |   \)"<<std::endl;
 std::cout <<R"(|____|__ \|___|_______ \_______ \/_______  / |____|_  /    |_______ \_______  /_______  /_______  //_______  / |____|_  /)"<<std::endl;
 std::cout <<R"(        \/            \/       \/        \/         \/             \/       \/        \/        \/         \/         \/)"<<std::endl;
 std::cout <<R"(=========================================================================================================================)"<<std::endl;
  // read in the first argument
  std::string filename = argv[1];
  
  // print message
  std::cout << "Filename: \"" << filename << "\"" << std::endl;
  Settings settings;
  settings.loadFromFile(filename);
  settings.printSettings();
 
  // write 5 output files
  for (int i = 0; i < 5; i++)
  {
    writeParaviewOutput(i);
  }

  return EXIT_SUCCESS;
}