
#include "computation/computation.h"
#include "computation/computation_parallel.h"
#include <iostream>
#include <cstdlib>
#include <mpi.h>

int main(int argc, char *argv[])
{
  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }
  // read in the first argument
  std::string filename = argv[1];
  MPI_Init(NULL, NULL);
  auto computation = ComputationParallel();
  computation.initialize(filename);
  computation.runSimulation();
  MPI_Finalize();
  return EXIT_SUCCESS;
}