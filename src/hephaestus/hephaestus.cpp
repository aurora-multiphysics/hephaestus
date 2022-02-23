//            -----------------------------------------------------
//            Hephaestus Miniapp: Electromagnetics
//            -----------------------------------------------------
//
// This miniapp is intended to provide a configurable interface to MFEM-based
// finite element solvers for electromagnetics problems.

#include <memory>
#include <iostream>
#include <fstream>
#include "hephaestus.hpp"

using namespace std;
using namespace mfem;
using namespace mfem::common;
using namespace mfem::electromagnetics;

double potential(const mfem::Vector &x, double t)
{
   double wj_ (2.0*M_PI/60.0);
   // the value
   double T;
   if (x[2] < 0.0)
   {
      T = 1.0;
   }
   else
   {
      T = -1.0;
   }

   return T*cos(wj_ * t);
}

hephaestus::Inputs joule_example_inputs()
{
   hephaestus::BCMap bc_map;
   bc_map.setBC(std::string("tangential_dEdt"), hephaestus::BoundaryCondition(std::string("boundary_1"), Array<int>({1,2,3})));
   bc_map.setBC(std::string("thermal_flux"), hephaestus::BoundaryCondition(std::string("boundary_2"), Array<int>({1,2})));
   hephaestus::BoundaryCondition poisson_bc(std::string("boundary_3"), Array<int>({1,2}));
   poisson_bc.scalar_func = potential;
   bc_map.setBC(std::string("electric_potential"), poisson_bc);

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

   hephaestus::Inputs inputs(std::string("cylinder-hex-q2.gen"), std::string("Joule"), 2, bc_map, material_map);
   return inputs;
}

int main(int argc, char *argv[])
{
   MPI_Session mpi(argc, argv);
   int myid;
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   // Parse command line arguments
   const char *formulation = "None";
   mfem::OptionsParser args(argc, argv);
   args.AddOption(&formulation, "-form", "--formulation",
                  "Name of formulation to use during solve.");
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(cout);
      }
      return 1;
   }

   // Create example inputs based on formulation name
   hephaestus::Inputs inputs;
   if (strcmp(formulation,"Joule")==0)
   {
      inputs = joule_example_inputs();
   }
   else if (strcmp(formulation,"None")==0)
   {
      std::cout<<"Formulation name not provided. \n";
      exit(0);
   }
   else
   {
      std::cout<<"Formulation name " << formulation << " not recognised. \n";
      exit(0);
   }

   // Launch
   run_hephaestus(argc, argv, inputs);

}