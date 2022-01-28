//            -----------------------------------------------------
//            Hephaestus Miniapp: Electromagnetics
//            -----------------------------------------------------
//
// This miniapp is intended to provide a configurable interface to MFEM-based
// finite element solvers for electromagnetics problems.

#include <memory>
#include <iostream>
#include <fstream>
#include "hephaestus_joule.hpp"

using namespace std;
using namespace mfem;
using namespace mfem::common;
using namespace mfem::electromagnetics;


int main(int argc, char *argv[])
{
   MPI_Session mpi(argc, argv);
   joule_solve(argc, argv);
}