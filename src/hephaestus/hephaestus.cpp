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

   std::vector<BCMap> bc_maps({
      BCMap(std::string("curl_bc"), Array<int>({1,2,3})),
      BCMap(std::string("thermal_bc"), Array<int>({1,2})),
      BCMap(std::string("poisson_bc"), Array<int>({1,2})),
   });

   Inputs inputs(bc_maps);
   joule_solve(argc, argv, inputs);
}