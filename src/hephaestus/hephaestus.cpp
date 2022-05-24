// //            -----------------------------------------------------
// //            Hephaestus Miniapp: Electromagnetics
// //            -----------------------------------------------------
// //
// // This miniapp is intended to provide a configurable interface to MFEM-based
// // finite element solvers for electromagnetics problems.

// #include "hephaestus.hpp"

// #include <fstream>
// #include <iostream>
// #include <memory>

// using namespace std;
// using namespace mfem;
// using namespace mfem::common;
// using namespace mfem::electromagnetics;

// double potential(const mfem::Vector &x, double t) {
//   double wj_(2.0 * M_PI / 60.0);
//   // the value
//   double T;
//   if (x[2] < 0.0) {
//     T = 1.0;
//   } else {
//     T = -1.0;
//   }

//   return T * cos(wj_ * t);
// }

// hephaestus::Inputs joule_example_inputs() {
//   hephaestus::BCMap bc_map;
//   bc_map["tangential_dEdt"] = new hephaestus::BoundaryCondition(
//       std::string("boundary_1"), Array<int>({1, 2, 3}));

//   bc_map["thermal_flux"] = new hephaestus::BoundaryCondition(
//       std::string("boundary_2"), Array<int>({1, 2}));

//   bc_map["electric_potential"] = new hephaestus::FunctionDirichletBC(
//       std::string("boundary_3"), Array<int>({1, 2}),
//       new mfem::FunctionCoefficient(potential));

//   double sigma = 2.0 * M_PI * 10;
//   double Tcapacity = 1.0;
//   double Tconductivity = 0.01;

//   double sigmaAir;
//   double TcondAir;
//   double TcapAir;

//   sigmaAir = 1.0e-6 * sigma;
//   TcondAir = 1.0e6 * Tconductivity;
//   TcapAir = 1.0 * Tcapacity;

//   hephaestus::Subdomain wire("wire", 1);
//   wire.property_map["electrical_conductivity"] = new ConstantCoefficient(sigma);
//   wire.property_map["heat_capacity"] = new ConstantCoefficient(Tcapacity);
//   wire.property_map["inverse_heat_capacity"] =
//       new ConstantCoefficient(1.0 / Tcapacity);
//   wire.property_map["inverse_thermal_conductivity"] =
//       new ConstantCoefficient(1.0 / Tconductivity);

//   hephaestus::Subdomain air("air", 2);
//   air.property_map["electrical_conductivity"] =
//       new ConstantCoefficient(sigmaAir);
//   air.property_map["heat_capacity"] = new ConstantCoefficient(TcapAir);
//   air.property_map["inverse_heat_capacity"] =
//       new ConstantCoefficient(1.0 / TcapAir);
//   air.property_map["inverse_thermal_conductivity"] =
//       new ConstantCoefficient(1.0 / TcondAir);

//   hephaestus::DomainProperties material_map(
//       std::vector<hephaestus::Subdomain>({wire, air}));

//   hephaestus::Executioner executioner(std::string("transient"), 0.5, 100.0);
//   mfem::Mesh mesh(std::string("cylinder-hex-q2.gen").c_str(), 1, 1);
//   hephaestus::Inputs inputs(mesh, std::string("Joule"), 2, bc_map, material_map,
//                             executioner);
//   return inputs;
// }

// void e_bc_r(const Vector &x, Vector &E) {
//   E.SetSize(3);
//   E = 0.0;
// }

// void e_bc_i(const Vector &x, Vector &E) {
//   E.SetSize(3);
//   E = 0.0;
// }

// hephaestus::Inputs hertz_example_inputs() {
//   hephaestus::BCMap bc_map;

//   // dirichlet
//   hephaestus::VectorFunctionDirichletBC e_bc(
//       std::string("boundary_1"), Array<int>({1, 2}),
//       new VectorFunctionCoefficient(3, e_bc_r),
//       new VectorFunctionCoefficient(3, e_bc_i));

//   bc_map["tangential_E"] = new hephaestus::VectorFunctionDirichletBC(e_bc);

//   hephaestus::Subdomain air("air", 1);

//   air.property_map["real_electrical_conductivity"] =
//       new ConstantCoefficient(0.0);
//   air.property_map["imag_electrical_conductivity"] =
//       new ConstantCoefficient(0.0);
//   air.property_map["real_rel_permittivity"] = new ConstantCoefficient(1.0);
//   air.property_map["imag_rel_permittivity"] = new ConstantCoefficient(0.0);
//   air.property_map["real_rel_permeability"] = new ConstantCoefficient(1.0);
//   air.property_map["imag_rel_permeability"] = new ConstantCoefficient(0.0);

//   hephaestus::DomainProperties material_map(
//       std::vector<hephaestus::Subdomain>({air}));

//   hephaestus::Executioner executioner(std::string("transient"), 0.5, 100.0);
//   mfem::Mesh mesh(std::string("irises.g").c_str(), 1, 1);
//   hephaestus::Inputs inputs(mesh, std::string("Hertz"), 2, bc_map, material_map,
//                             executioner);
//   return inputs;
// }

int main(int argc, char *argv[]) {
  // MPI_Session mpi(argc, argv);
  // int myid;
  // MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // // Parse command line arguments
  // const char *formulation = "None";
  // mfem::OptionsParser args(argc, argv);
  // args.AddOption(&formulation, "-form", "--formulation",
  //                "Name of formulation to use during solve.");
  // args.Parse();
  // if (!args.Good()) {
  //   if (myid == 0) {
  //     args.PrintUsage(cout);
  //   }
  //   return 1;
  // }

  // // Create example inputs based on formulation name
  // hephaestus::Inputs inputs;
  // if (strcmp(formulation, "Joule") == 0) {
  //   inputs = joule_example_inputs();
  // } else if (strcmp(formulation, "Hertz") == 0) {
  //   inputs = hertz_example_inputs();
  // } else if (strcmp(formulation, "None") == 0) {
  //   std::cout << "Formulation name not provided. \n";
  //   exit(0);
  // } else {
  //   std::cout << "Formulation name " << formulation << " not recognised. \n";
  //   exit(0);
  // }

  // // Launch
  // run_hephaestus(argc, argv, inputs);
}