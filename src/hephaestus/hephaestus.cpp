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

   hephaestus::BCMap bc_map;
   bc_map.setBC(std::string("curl_bc"), hephaestus::BoundaryCondition(std::string("boundary_1"), Array<int>({1,2,3})));
   bc_map.setBC(std::string("thermal_bc"), hephaestus::BoundaryCondition(std::string("boundary_2"), Array<int>({1,2})));
   bc_map.setBC(std::string("poisson_bc"), hephaestus::BoundaryCondition(std::string("boundary_3"), Array<int>({1,2})));

   // std::vector<hephaestus::BoundaryCondition> bc_maps({
   //    hephaestus::BoundaryCondition(std::string("curl_bc"), Array<int>({1,2,3})),
   //    hephaestus::BoundaryCondition(std::string("thermal_bc"), Array<int>({1,2})),
   //    hephaestus::BoundaryCondition(std::string("poisson_bc"), Array<int>({1,2})),
   // });

   double sigma = 2.0*M_PI*10;
   double Tcapacity = 1.0;
   double Tconductivity = 0.01;

   double sigmaAir;
   double TcondAir;
   double TcapAir;

   sigmaAir     = 1.0e-6 * sigma;
   TcondAir     = 1.0e6  * Tconductivity;
   TcapAir      = 1.0    * Tcapacity;

   hephaestus::Material copper("copper", 1);
   copper.setMaterialProperty(std::string("electrical_conductivity"), sigma);
   copper.setMaterialProperty(std::string("thermal_conductivity"), Tconductivity);
   copper.setMaterialProperty(std::string("heat_capacity"), Tcapacity);

   hephaestus::Material air("air", 2);
   air.setMaterialProperty(std::string("electrical_conductivity"), sigmaAir);
   air.setMaterialProperty(std::string("thermal_conductivity"), TcondAir);
   air.setMaterialProperty(std::string("heat_capacity"), TcapAir);

   hephaestus::MaterialMap material_map(std::vector<hephaestus::Material>({copper, air}));

   hephaestus::Inputs inputs(std::string("rod"), std::string("cylinder-hex-q2.gen"), bc_map, material_map);
   joule_solve(argc, argv, inputs);
}