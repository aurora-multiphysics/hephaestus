#include "hephaestus.hpp"

const char *DATA_DIR = "../../data/";

hephaestus::Coefficients defineCoefficients() {
  hephaestus::Subdomain air("air", 1);
  air.scalar_coefficients.Register("electrical_conductivity",
                                   new mfem::ConstantCoefficient(1.0), true);
  hephaestus::Subdomain plate("plate", 2);
  plate.scalar_coefficients.Register(
      "electrical_conductivity", new mfem::ConstantCoefficient(3.526e7), true);
  hephaestus::Subdomain coil1("coil1", 3);
  coil1.scalar_coefficients.Register("electrical_conductivity",
                                     new mfem::ConstantCoefficient(1.0), true);
  hephaestus::Subdomain coil2("coil2", 4);
  coil2.scalar_coefficients.Register("electrical_conductivity",
                                     new mfem::ConstantCoefficient(1.0), true);
  hephaestus::Subdomain coil3("coil3", 5);
  coil3.scalar_coefficients.Register("electrical_conductivity",
                                     new mfem::ConstantCoefficient(1.0), true);
  hephaestus::Subdomain coil4("coil4", 6);
  coil4.scalar_coefficients.Register("electrical_conductivity",
                                     new mfem::ConstantCoefficient(1.0), true);
  hephaestus::Coefficients coefficients(std::vector<hephaestus::Subdomain>(
      {air, plate, coil1, coil2, coil3, coil4}));

  coefficients.scalars.Register("frequency",
                                new mfem::ConstantCoefficient(200.0), true);
  coefficients.scalars.Register("dielectric_permittivity",
                                new mfem::ConstantCoefficient(8.854e-12), true);
  coefficients.scalars.Register("magnetic_permeability",
                                new mfem::ConstantCoefficient(M_PI * 4.0e-7),
                                true);

  double Itotal = 2742;
  coefficients.scalars.Register("I", new mfem::ConstantCoefficient(Itotal),
                                true);
  return coefficients;
}

hephaestus::Sources defineSources() {
  hephaestus::InputParameters div_free_source_params;
  // This vector of subdomains will form the coil that we pass to
  // ClosedCoilSolver
  int order = 1;
  int electrode_attr = 7;
  std::string coil_attr = "3 4 5 6";
  mfem::Array<int> coil_domains;
  std::stringstream ss(coil_attr);
  int att;
  while (ss >> att)
    coil_domains.Append(att);

  hephaestus::InputParameters coilsolver_pars;
  coilsolver_pars.SetParam("HCurlFESpaceName", std::string("HCurl"));
  coilsolver_pars.SetParam("JGridFunctionName",
                           std::string("source_current_density"));
  coilsolver_pars.SetParam("IFuncCoefName", std::string("I"));
  coilsolver_pars.SetParam("H1FESpaceName", std::string("H1"));
  coilsolver_pars.SetParam("JTransfer", true);

  hephaestus::Sources sources;
  sources.Register("source",
                   new hephaestus::ClosedCoilSolver(
                       coilsolver_pars, coil_domains, electrode_attr),
                   true);
  return sources;
}

hephaestus::Outputs defineOutputs() {
  hephaestus::Outputs outputs;
  outputs.Register("ParaViewDataCollection",
                   new mfem::ParaViewDataCollection("ComplexTeam7ParaView"),
                   true);
  return outputs;
}

int main(int argc, char *argv[]) {
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&DATA_DIR, "-dataDir", "--data_directory",
                 "Directory storing input data for tests.");
  args.Parse();
  MPI_Init(&argc, &argv);

  // Create Formulation
  hephaestus::ComplexAFormulation *problem_builder =
      new hephaestus::ComplexAFormulation(
          "magnetic_reluctivity", "electrical_conductivity",
          "dielectric_permittivity", "frequency", "magnetic_vector_potential",
          "magnetic_vector_potential_real", "magnetic_vector_potential_imag");
  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./team7.g")).c_str(), 1,
                  1);
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(mfem::ParMesh(MPI_COMM_WORLD, mesh));
  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace("HCurl", "ND_3D_P1");
  problem_builder->AddFESpace("H1", "H1_3D_P1");
  problem_builder->AddFESpace("HDiv", "RT_3D_P0");
  problem_builder->AddFESpace("L2", "L2_3D_P0");

  problem_builder->AddGridFunction("magnetic_vector_potential_real", "HCurl");
  problem_builder->AddGridFunction("magnetic_vector_potential_imag", "HCurl");

  problem_builder->AddGridFunction("electric_field_real", "HCurl");
  problem_builder->AddGridFunction("electric_field_imag", "HCurl");
  problem_builder->registerElectricFieldAux("electric_field_real",
                                            "electric_field_imag");

  problem_builder->AddGridFunction("magnetic_flux_density_real", "HDiv");
  problem_builder->AddGridFunction("magnetic_flux_density_imag", "HDiv");
  problem_builder->registerMagneticFluxDensityAux("magnetic_flux_density_real",
                                                  "magnetic_flux_density_imag");

  problem_builder->AddGridFunction("current_density_real", "HDiv");
  problem_builder->AddGridFunction("current_density_imag", "HDiv");
  problem_builder->registerCurrentDensityAux("current_density_real",
                                             "current_density_imag");

  problem_builder->AddGridFunction("joule_heating_density", "L2");
  problem_builder->registerJouleHeatingDensityAux(
      "joule_heating_density", "electric_field_real", "electric_field_imag",
      "electrical_conductivity");

  hephaestus::Coefficients coefficients = defineCoefficients();
  problem_builder->SetCoefficients(coefficients);

  hephaestus::Sources sources = defineSources();
  problem_builder->SetSources(sources);

  hephaestus::Outputs outputs = defineOutputs();
  problem_builder->SetOutputs(outputs);

  hephaestus::InputParameters solver_options;
  solver_options.SetParam("Tolerance", float(1.0e-16));
  solver_options.SetParam("MaxIter", (unsigned int)1000);
  solver_options.SetParam("PrintLevel", 0);
  problem_builder->SetSolverOptions(solver_options);

  hephaestus::ProblemBuildSequencer sequencer(problem_builder);
  sequencer.ConstructOperatorProblem();
  std::unique_ptr<hephaestus::SteadyStateProblem> problem =
      problem_builder->ReturnProblem();
  problem.get()->outputs.EnableGLVis(true);

  hephaestus::InputParameters exec_params;
  exec_params.SetParam("Problem", problem.get());
  hephaestus::SteadyExecutioner *executioner =
      new hephaestus::SteadyExecutioner(exec_params);

  std::cout << "Created exec ";
  executioner->Execute();

  MPI_Finalize();
}
