#include "hephaestus.hpp"

const char *DATA_DIR = "../../data/";

static void source_current(const mfem::Vector &xv, double t, mfem::Vector &J) {
  double x0(194e-3);  // Coil centre x coordinate
  double y0(100e-3);  // Coil centre y coordinate
  double a(50e-3);    // Coil thickness
  double I0(2742);    // Coil current in Ampere-turns
  double S(2.5e-3);   // Coil cross sectional area
  double freq(200.0); // Frequency in Hz

  double x = xv(0);
  double y = xv(1);

  // Current density magnitude
  double Jmag = (I0 / S);

  // Calculate x component of current density unit vector
  if (abs(x - x0) < a) {
    J(0) = -(y - y0) / abs(y - y0);
  } else if (abs(y - y0) < a) {
    J(0) = 0.0;
  } else {
    J(0) = -(y - (y0 + a * ((y - y0) / abs(y - y0)))) /
           hypot(x - (x0 + a * ((x - x0) / abs(x - x0))),
                 y - (y0 + a * ((y - y0) / abs(y - y0))));
  }

  // Calculate y component of current density unit vector
  if (abs(y - y0) < a) {
    J(1) = (x - x0) / abs(x - x0);
  } else if (abs(x - x0) < a) {
    J(1) = 0.0;
  } else {
    J(1) = (x - (x0 + a * ((x - x0) / abs(x - x0)))) /
           hypot(x - (x0 + a * ((x - x0) / abs(x - x0))),
                 y - (y0 + a * ((y - y0) / abs(y - y0))));
  }

  // Calculate z component of current density unit vector
  J(2) = 0.0;

  // Scale by current density magnitude
  J *= Jmag;
}

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
  coefficients.scalars.Register("magnetic_permeability",
                                new mfem::ConstantCoefficient(M_PI * 4.0e-7),
                                true);

  mfem::VectorFunctionCoefficient *JSrcCoef =
      new mfem::VectorFunctionCoefficient(3, source_current);
  mfem::Array<mfem::VectorCoefficient *> sourcecoefs(4);
  sourcecoefs[0] = JSrcCoef;
  sourcecoefs[1] = JSrcCoef;
  sourcecoefs[2] = JSrcCoef;
  sourcecoefs[3] = JSrcCoef;
  mfem::Array<int> coilsegments(4);
  coilsegments[0] = 3;
  coilsegments[1] = 4;
  coilsegments[2] = 5;
  coilsegments[3] = 6;
  mfem::PWVectorCoefficient *JSrcRestricted =
      new mfem::PWVectorCoefficient(3, coilsegments, sourcecoefs);
  coefficients.vectors.Register("source", JSrcRestricted, true);
  return coefficients;
}

hephaestus::Sources defineSources() {
  hephaestus::InputParameters div_free_source_params;
  div_free_source_params.SetParam("SourceName", std::string("source"));
  div_free_source_params.SetParam("HCurlFESpaceName", std::string("HCurl"));
  div_free_source_params.SetParam("H1FESpaceName", std::string("H1"));
  hephaestus::InputParameters current_solver_options;
  current_solver_options.SetParam("Tolerance", float(1.0e-12));
  current_solver_options.SetParam("MaxIter", (unsigned int)200);
  current_solver_options.SetParam("PrintLevel", 0);
  div_free_source_params.SetParam("SolverOptions", current_solver_options);
  hephaestus::Sources sources;
  sources.Register("source",
                   new hephaestus::DivFreeSource(div_free_source_params), true);
  return sources;
}
hephaestus::Outputs defineOutputs() {
  std::map<std::string, mfem::DataCollection *> data_collections;
  data_collections["ParaViewDataCollection"] =
      new mfem::ParaViewDataCollection("ComplexTeam7ParaView");
  hephaestus::Outputs outputs(data_collections);
  return outputs;
}

int main(int argc, char *argv[]) {
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&DATA_DIR, "-dataDir", "--data_directory",
                 "Directory storing input data for tests.");
  args.Parse();
  MPI_Init(&argc, &argv);

  // Create Formulation
  hephaestus::SteadyStateProblemBuilder *problem_builder =
      new hephaestus::ComplexAFormulation(
          "magnetic_reluctivity", "electrical_conductivity",
          "dielectric_permittivity", "frequency", "magnetic_vector_potential");
  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./team7.g")).c_str(), 1,
                  1);
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(mfem::ParMesh(MPI_COMM_WORLD, mesh));
  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P1"));
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P1"));
  problem_builder->AddFESpace(std::string("HDiv"), std::string("RT_3D_P0"));
  problem_builder->AddFESpace(std::string("L2"), std::string("L2_3D_P1"));
  problem_builder->AddGridFunction(
      std::string("magnetic_vector_potential_real"), std::string("HCurl"));
  problem_builder->AddGridFunction(
      std::string("magnetic_vector_potential_imag"), std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("electric_field_real"),
                                   std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("electric_field_imag"),
                                   std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("magnetic_flux_density_real"),
                                   std::string("HDiv"));
  problem_builder->AddGridFunction(std::string("magnetic_flux_density_imag"),
                                   std::string("HDiv"));
  problem_builder->AddGridFunction(std::string("joule_heating_load"),
                                   std::string("L2"));
  hephaestus::Coefficients coefficients = defineCoefficients();
  problem_builder->SetCoefficients(coefficients);

  hephaestus::Sources sources = defineSources();
  problem_builder->SetSources(sources);

  hephaestus::Outputs outputs = defineOutputs();
  problem_builder->SetOutputs(outputs);

  hephaestus::AuxSolvers postprocessors;
  hephaestus::InputParameters joule_heating_params;
  joule_heating_params.SetParam("VariableName",
                                std::string("joule_heating_load"));
  joule_heating_params.SetParam("CoefficientName", std::string("JouleHeating"));
  joule_heating_params.SetParam("ElectricFieldName",
                                std::string("magnetic_vector_potential"));
  joule_heating_params.SetParam("ConductivityCoefName",
                                std::string("electrical_conductivity"));
  joule_heating_params.SetParam("ComplexField", true);
  postprocessors.Register(
      "JouleHeatingAux",
      new hephaestus::JouleHeatingAuxSolver(joule_heating_params), true);
  problem_builder->SetPostprocessors(postprocessors);

  hephaestus::InputParameters solver_options;
  solver_options.SetParam("Tolerance", float(1.0e-16));
  solver_options.SetParam("MaxIter", (unsigned int)1000);
  solver_options.SetParam("PrintLevel", 0);
  problem_builder->SetSolverOptions(solver_options);

  hephaestus::ProblemBuildSequencer sequencer(problem_builder);
  sequencer.ConstructOperatorProblem();
  std::unique_ptr<hephaestus::SteadyStateProblem> problem =
      problem_builder->ReturnProblem();

  hephaestus::InputParameters exec_params;
  exec_params.SetParam("UseGLVis", true);
  exec_params.SetParam("Problem", problem.get());
  hephaestus::SteadyExecutioner *executioner =
      new hephaestus::SteadyExecutioner(exec_params);

  std::cout << "Created exec ";
  executioner->Init();
  executioner->Execute();

  MPI_Finalize();
}
