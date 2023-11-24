#include "hephaestus.hpp"

const char *DATA_DIR = "../../data/";

hephaestus::Coefficients defineCoefficients(double Itotal) {

  hephaestus::Subdomain coil("coil", 1);
  coil.scalar_coefficients.Register(
      "electrical_conductivity", new mfem::ConstantCoefficient(3.526e7), true);
  hephaestus::Subdomain air("air", 2);
  air.scalar_coefficients.Register("electrical_conductivity",
                                   new mfem::ConstantCoefficient(1.0), true);
  hephaestus::Coefficients coefficients(
      std::vector<hephaestus::Subdomain>({coil, air}));
  coefficients.scalars.Register("magnetic_permeability",
                                new mfem::ConstantCoefficient(M_PI * 4.0e-7),
                                true);

  // Time-dependent current
  coefficients.scalars.Register("I", new mfem::ConstantCoefficient(Itotal),
                                true);

  return coefficients;
}

hephaestus::Sources defineSources(std::pair<int, int> elec,
                                  mfem::Array<int> coil_domains) {

  hephaestus::InputParameters coilsolver_pars;
  coilsolver_pars.SetParam("SourceName", std::string("source_current_density"));
  coilsolver_pars.SetParam("PotentialName", std::string("auxiliary_potential"));
  coilsolver_pars.SetParam("IFuncCoefName", std::string("I"));
  coilsolver_pars.SetParam("HelmholtzProjection", false);

  hephaestus::Sources sources;
  sources.Register(
      "source",
      new hephaestus::OpenCoilSolver(coilsolver_pars, coil_domains, elec),
      true);
  return sources;
}

hephaestus::Outputs defineOutputs() {
  std::map<std::string, mfem::DataCollection *> data_collections;
  data_collections["ParaViewDataCollection"] =
      new mfem::ParaViewDataCollection("Team7ParaView");
  hephaestus::Outputs outputs(data_collections);
  return outputs;
}

int main(int argc, char *argv[]) {

  // Refinement and order
  int par_ref_lvl = -1;
  int order = 1;

  // Total electrical current going around the coil. Must be nonzero, can be
  // changed later.
  double Itotal = 10;

  // Attribute that defines the internal faces over which we apply the potential
  // difference
  std::pair<int, int> elec_attrs;

  // Mesh file
  std::string mesh_filename = "coil.gen";

  // Domain attributes of the coil to be solved and boundary attributes of the
  // electrodes
  std::string coil_attr = "1";
  std::string elec_bdr_attr = "1 2";

  mfem::OptionsParser args(argc, argv);
  args.AddOption(&DATA_DIR, "-dataDir", "--data_directory",
                 "Directory storing input data for tests.");
  args.AddOption(&par_ref_lvl, "-ref", "--parallel-refinement",
                 "Parallel refinement level.");
  args.AddOption(&order, "-o", "--order", "Base functions order");
  args.AddOption(&Itotal, "-I", "--Itotal", "Total electrical current.");
  args.AddOption(&mesh_filename, "-f", "--mesh-filename", "Mesh file name");
  args.AddOption(
      &coil_attr, "-cd", "--coil-domains",
      "List of coil domain attributes separated by spaces, e.g. \'1 3 4\'");
  args.AddOption(&elec_bdr_attr, "-e", "--electrode-attrs",
                 "List of electrode attributes separated by spaces, e.g. \'1 "
                 "2\'. Must be two values.");

  // ADD AN OPTION FOR DEFINING THE BOUNDARY ATTRIBUTES OF THE ELECTRODES ////
  args.Parse();

  MPI_Init(&argc, &argv);

  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + mesh_filename).c_str(), 1, 1);
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(mfem::ParMesh(MPI_COMM_WORLD, mesh));
  for (int l = 0; l < par_ref_lvl; ++l)
    pmesh->UniformRefinement();

  // This vector of subdomains will form the coil that we pass to
  // OpenCoilSolver
  mfem::Array<int> coil_domains;

  // This vector of attributes will form the electrodes that we pass to
  // OpenCoilSolver
  mfem::Array<int> elec_bdr_array;

  // Parsing the string of attributes
  std::stringstream ss(coil_attr);
  int att;
  while (ss >> att)
    coil_domains.Append(att);

  std::stringstream ss_bdr(elec_bdr_attr);
  while (ss_bdr >> att)
    elec_bdr_array.Append(att);

  if (elec_bdr_array.Size() != 2)
    mfem::mfem_error(
        "Electrode boundary attribute list must contain two attributes.");

  elec_attrs.first = elec_bdr_array[0];
  elec_attrs.second = elec_bdr_array[1];

  // Create Formulation
  hephaestus::MagnetostaticFormulation *problem_builder =
      new hephaestus::MagnetostaticFormulation("magnetic_reluctivity",
                                               "magnetic_permeability",
                                               "magnetic_vector_potential");
  // Set Mesh
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
  hephaestus::Coefficients coefficients = defineCoefficients(Itotal);
  problem_builder->SetCoefficients(coefficients);

  hephaestus::Sources sources = defineSources(elec_attrs, coil_domains);
  problem_builder->SetSources(sources);

  hephaestus::Outputs outputs = defineOutputs();
  problem_builder->SetOutputs(outputs);

  hephaestus::InputParameters solver_options;
  solver_options.SetParam("Tolerance", float(1.0e-18));
  solver_options.SetParam("AbsTolerance", float(1.0e-18));
  solver_options.SetParam("MaxIter", (unsigned int)1000);
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
  executioner->Init();
  executioner->Execute();

  MPI_Finalize();
}