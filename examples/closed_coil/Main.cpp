#include "hephaestus.hpp"

const char *DATA_DIR = "../../data/";

hephaestus::Coefficients defineCoefficients() {
  hephaestus::Coefficients coefficients;
  coefficients.scalars.Register("magnetic_permeability",
                                new mfem::ConstantCoefficient(M_PI * 4.0e-7),
                                true);
  coefficients.scalars.Register("electrical_conductivity",
                                new mfem::ConstantCoefficient(1.0), true);

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
                   new mfem::ParaViewDataCollection("ClosedCoilParaView"),
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
  hephaestus::MagnetostaticFormulation *problem_builder =
      new hephaestus::MagnetostaticFormulation("magnetic_reluctivity",
                                               "magnetic_permeability",
                                               "magnetic_vector_potential");
  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./team7.g")).c_str(), 1,
                  1);
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(mfem::ParMesh(MPI_COMM_WORLD, mesh));

  int par_ref_lvl = 1;
  for (int l = 0; l < par_ref_lvl; ++l)
    pmesh->UniformRefinement();

  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P1"));
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P1"));
  problem_builder->AddFESpace(std::string("HDiv"), std::string("RT_3D_P0"));
  problem_builder->AddGridFunction(std::string("magnetic_vector_potential"),
                                   std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("source_current_density"),
                                   std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("magnetic_flux_density"),
                                   std::string("HDiv"));
  problem_builder->registerMagneticFluxDensityAux("magnetic_flux_density");
  hephaestus::Coefficients coefficients = defineCoefficients();
  problem_builder->SetCoefficients(coefficients);

  hephaestus::Sources sources = defineSources();
  problem_builder->SetSources(sources);

  hephaestus::Outputs outputs = defineOutputs();
  problem_builder->SetOutputs(outputs);

  hephaestus::InputParameters solver_options;
  solver_options.SetParam("Tolerance", float(1.0e-13));
  solver_options.SetParam("AbsTolerance", float(1.0e-16));
  solver_options.SetParam("MaxIter", (unsigned int)500);
  solver_options.SetParam("PrintLevel", 2);
  problem_builder->SetSolverOptions(solver_options);

  hephaestus::ProblemBuildSequencer sequencer(problem_builder);
  sequencer.ConstructEquationSystemProblem();
  std::unique_ptr<hephaestus::SteadyStateProblem> problem =
      problem_builder->ReturnProblem();
  hephaestus::InputParameters exec_params;
  exec_params.SetParam("VisualisationSteps", int(1));
  exec_params.SetParam("UseGLVis", true);
  exec_params.SetParam("Problem", problem.get());
  hephaestus::SteadyExecutioner *executioner =
      new hephaestus::SteadyExecutioner(exec_params);

  mfem::out << "Created executioner";
  executioner->Execute();

  MPI_Finalize();
}