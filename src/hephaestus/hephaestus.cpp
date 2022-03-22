//            -----------------------------------------------------
//            Hephaestus Miniapp: Electromagnetics
//            -----------------------------------------------------
//
// This miniapp is intended to provide a configurable interface to MFEM-based
// finite element solvers for electromagnetics problems.

#include "hephaestus.hpp"

#include <fstream>
#include <iostream>
#include <memory>

using namespace std;
using namespace mfem;
using namespace mfem::common;
using namespace mfem::electromagnetics;

double potential(const mfem::Vector &x, double t) {
  double wj_(2.0 * M_PI / 60.0);
  // the value
  double T;
  if (x[2] < 0.0) {
    T = 1.0;
  } else {
    T = -1.0;
  }

  return T * cos(wj_ * t);
}

hephaestus::Inputs joule_example_inputs() {
  hephaestus::BCMap bc_map;
  bc_map["tangential_dEdt"] = new hephaestus::BoundaryCondition(
      std::string("boundary_1"), Array<int>({1, 2, 3}));

  bc_map["thermal_flux"] = new hephaestus::BoundaryCondition(
      std::string("boundary_2"), Array<int>({1, 2}));

  bc_map["electric_potential"] = new hephaestus::FunctionDirichletBC(
      std::string("boundary_3"), Array<int>({1, 2}),
      new mfem::FunctionCoefficient(potential));

  double sigma = 2.0 * M_PI * 10;
  double Tcapacity = 1.0;
  double Tconductivity = 0.01;

  double sigmaAir;
  double TcondAir;
  double TcapAir;

  sigmaAir = 1.0e-6 * sigma;
  TcondAir = 1.0e6 * Tconductivity;
  TcapAir = 1.0 * Tcapacity;

  hephaestus::Material copper("copper", 1);
  copper.setMaterialProperty(std::string("electrical_conductivity"), sigma);
  copper.setMaterialProperty(std::string("thermal_conductivity"),
                             Tconductivity);
  copper.setMaterialProperty(std::string("heat_capacity"), Tcapacity);

  hephaestus::Material air("air", 2);
  air.setMaterialProperty(std::string("electrical_conductivity"), sigmaAir);
  air.setMaterialProperty(std::string("thermal_conductivity"), TcondAir);
  air.setMaterialProperty(std::string("heat_capacity"), TcapAir);

  hephaestus::MaterialMap material_map(
      std::vector<hephaestus::Material>({copper, air}));

  hephaestus::Executioner executioner(std::string("transient"), 0.5, 100.0);
  hephaestus::Inputs inputs(std::string("cylinder-hex-q2.gen"),
                            std::string("Joule"), 2, bc_map, material_map,
                            executioner);
  return inputs;
}

void e_bc_r(const Vector &x, Vector &E) {
  E.SetSize(3);
  E = 0.0;
}

void e_bc_i(const Vector &x, Vector &E) {
  E.SetSize(3);
  E = 0.0;
}

hephaestus::Inputs hertz_example_inputs() {
  hephaestus::BCMap bc_map;

  // dirichlet
  hephaestus::VectorFunctionDirichletBC e_bc(std::string("boundary_1"),
                                             Array<int>({1, 2}));
  e_bc.vector_func = e_bc_r;
  e_bc.vector_func_im = e_bc_i;

  bc_map["tangential_E"] = new hephaestus::VectorFunctionDirichletBC(e_bc);

  hephaestus::Material air("air", 1);
  air.setMaterialProperty(std::string("real_electrical_conductivity"), 0.0);
  air.setMaterialProperty(std::string("imag_electrical_conductivity"), 0.0);
  air.setMaterialProperty(std::string("real_rel_permittivity"), 1.0);
  air.setMaterialProperty(std::string("imag_rel_permittivity"), 0.0);
  air.setMaterialProperty(std::string("real_rel_permeability"), 1.0);
  air.setMaterialProperty(std::string("imag_rel_permeability"), 0.0);

  hephaestus::MaterialMap material_map(
      std::vector<hephaestus::Material>({air}));

  hephaestus::Executioner executioner(std::string("transient"), 0.5, 100.0);
  hephaestus::Inputs inputs(std::string("irises.g"), std::string("Hertz"), 2,
                            bc_map, material_map, executioner);
  return inputs;
}

int main(int argc, char *argv[]) {
  MPI_Session mpi(argc, argv);
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // Parse command line arguments
  const char *formulation = "None";
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&formulation, "-form", "--formulation",
                 "Name of formulation to use during solve.");
  args.Parse();
  if (!args.Good()) {
    if (myid == 0) {
      args.PrintUsage(cout);
    }
    return 1;
  }

  // Create example inputs based on formulation name
  hephaestus::Inputs inputs;
  if (strcmp(formulation, "Joule") == 0) {
    inputs = joule_example_inputs();
  } else if (strcmp(formulation, "Hertz") == 0) {
    inputs = hertz_example_inputs();
  } else if (strcmp(formulation, "None") == 0) {
    std::cout << "Formulation name not provided. \n";
    exit(0);
  } else {
    std::cout << "Formulation name " << formulation << " not recognised. \n";
    exit(0);
  }

  // Launch
  run_hephaestus(argc, argv, inputs);
}