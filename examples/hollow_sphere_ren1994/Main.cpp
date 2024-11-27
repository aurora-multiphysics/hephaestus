#include "hephaestus.hpp"

const char * DATA_DIR = "../data/";
const char * MESH_NAME = "hollow_sphere_vac.e";

hephaestus::Coefficients
defineCoefficients()
{
  double air_permeability = M_PI * 4.0e-7;
  double solid_permeability = air_permeability * 500.0;

  double solid_conductivity = 1.0;
  double air_conductivity = 1.0;

  hephaestus::Subdomain vacuum_region("vacuum_region", 107);
  vacuum_region._scalar_coefficients.Register(
      "electrical_conductivity", std::make_shared<mfem::ConstantCoefficient>(air_conductivity));
  vacuum_region._scalar_coefficients.Register(
      "magnetic_permeability", std::make_shared<mfem::ConstantCoefficient>(air_permeability));

  hephaestus::Subdomain sphere("sphere", 100);
  sphere._scalar_coefficients.Register(
      "electrical_conductivity", std::make_shared<mfem::ConstantCoefficient>(solid_conductivity));
  sphere._scalar_coefficients.Register(
      "magnetic_permeability", std::make_shared<mfem::ConstantCoefficient>(solid_permeability));

  hephaestus::Subdomain coil_0("coil_0", 103);
  coil_0._scalar_coefficients.Register(
      "electrical_conductivity", std::make_shared<mfem::ConstantCoefficient>(air_conductivity));
  coil_0._scalar_coefficients.Register(
      "magnetic_permeability", std::make_shared<mfem::ConstantCoefficient>(air_permeability));

  hephaestus::Subdomain coil_1("coil_1", 104);
  coil_1._scalar_coefficients.Register(
      "electrical_conductivity", std::make_shared<mfem::ConstantCoefficient>(air_conductivity));
  coil_1._scalar_coefficients.Register(
      "magnetic_permeability", std::make_shared<mfem::ConstantCoefficient>(air_permeability));

  hephaestus::Subdomain coil_2("coil_2", 105);
  coil_2._scalar_coefficients.Register(
      "electrical_conductivity", std::make_shared<mfem::ConstantCoefficient>(air_conductivity));
  coil_2._scalar_coefficients.Register(
      "magnetic_permeability", std::make_shared<mfem::ConstantCoefficient>(air_permeability));

  hephaestus::Subdomain coil_3("coil_3", 106);
  coil_3._scalar_coefficients.Register(
      "electrical_conductivity", std::make_shared<mfem::ConstantCoefficient>(air_conductivity));
  coil_3._scalar_coefficients.Register(
      "magnetic_permeability", std::make_shared<mfem::ConstantCoefficient>(air_permeability));

  hephaestus::Coefficients coefficients(
      std::vector<hephaestus::Subdomain>({vacuum_region, sphere, coil_0, coil_1, coil_2, coil_3}));
  // coefficients._scalars.Register("frequency",
  // std::make_shared<mfem::ConstantCoefficient>(200.0));
  // coefficients._scalars.Register("dielectric_permittivity",
  //                                std::make_shared<mfem::ConstantCoefficient>(8.854e-12));

  coefficients._scalars.Register("I", std::make_shared<mfem::ConstantCoefficient>(20000));

  return coefficients;
}

hephaestus::Sources
defineSources()
{
  hephaestus::InputParameters div_free_source_params;

  // This vector of subdomains will form the coil that we pass to
  // ClosedCoilSolver
  int order = 1;
  int electrode_attr = 110; // coil face
  std::string coil_attr = "103 104 105 106";
  mfem::Array<int> coil_domains;
  std::stringstream ss(coil_attr);
  int att;

  while (ss >> att)
  {
    coil_domains.Append(att);
    std::cout << "coil domains appending " << att << std::endl;
  }

  hephaestus::Sources sources;
  sources.Register("source",
                   std::make_shared<hephaestus::ClosedCoilSolver>("source_grad_phi",
                                                                  "HCurl",
                                                                  "H1",
                                                                  "I",
                                                                  "electrical_conductivity",
                                                                  coil_domains,
                                                                  electrode_attr,
                                                                  true));
  return sources;
}

hephaestus::Outputs
defineOutputs()
{
  hephaestus::Outputs outputs;
  outputs.Register("ParaViewDataCollection",
                   std::make_shared<mfem::ParaViewDataCollection>("HollowSphereParaView"));
  return outputs;
}

int
main(int argc, char * argv[])
{
  int ref_level(0);
  mfem::OptionsParser args(argc, argv);
  args.AddOption(
      &DATA_DIR, "-dataDir", "--data_directory", "Directory storing input data for tests.");
  args.AddOption(&MESH_NAME, "-meshName", "--mesh_name", "File name.");
  args.AddOption(
      &ref_level, "-ref", "--refinement_level", "number of uniform refinement iterations.");
  args.Parse();
  MPI_Init(&argc, &argv);

  hephaestus::logger.set_level(spdlog::level::info);

  // Create Formulation
  auto problem_builder = std::make_unique<hephaestus::MagnetostaticFormulation>(
      "magnetic_reluctivity", "magnetic_permeability", "magnetic_vector_potential");

  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + std::string(MESH_NAME)).c_str(), 1, 1);
  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);

  for (int l = 0; l < ref_level; ++l)
    pmesh->UniformRefinement();

  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P1"));
  problem_builder->AddFESpace(std::string("Vector_H1"), std::string("H1_3D_P1"), 3);
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P1"));
  problem_builder->AddFESpace(std::string("HDiv"), std::string("RT_3D_P0"));
  problem_builder->AddFESpace(std::string("Scalar_L2"), std::string("L2_3D_P0"));
  problem_builder->AddFESpace(std::string("Vector_L2"), std::string("L2_3D_P0"), 3);
  problem_builder->AddGridFunction(std::string("magnetic_vector_potential"), std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("source_grad_phi"), std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("magnetic_flux_density"), std::string("HDiv"));

  problem_builder->AddGridFunction(std::string("magnetic_force"), std::string("Vector_L2"));
  problem_builder->RegisterMagneticFluxDensityAux("magnetic_flux_density");

  int outer_sphere_id = 101;
  int inner_sphere_id = 102;
  int sphere_volume_id = 100;
  mfem::Array<int> boundary_marker;
  boundary_marker.Append(outer_sphere_id);
  boundary_marker.Append(inner_sphere_id);
  mfem::Array<int> volume_marker;
  volume_marker.Append(sphere_volume_id);

  hephaestus::Coefficients coefficients = defineCoefficients();
  problem_builder->SetCoefficients(coefficients);

  hephaestus::Sources sources = defineSources();
  problem_builder->SetSources(sources);

  hephaestus::Outputs outputs = defineOutputs();
  problem_builder->SetOutputs(outputs);

  auto magnetic_force_monitor = std::make_shared<hephaestus::MagneticForceAux>(
      "magnetic_flux_density", "magnetic_vector_potential", boundary_marker, "magnetic_force"
      // "magnetic_flux_density", "magnetic_vector_potential", volume_marker, "magnetic_force",
      // false
  );
  magnetic_force_monitor->SetPriority(2);
  problem_builder->AddPostprocessor("MaxwellStressMonitor", magnetic_force_monitor);

  hephaestus::InputParameters solver_options;
  solver_options.SetParam("Tolerance", float(1.0e-13));
  solver_options.SetParam("AbsTolerance", float(1.0e-16));
  solver_options.SetParam("MaxIter", (unsigned int)100);
  problem_builder->SetSolverOptions(solver_options);

  problem_builder->FinalizeProblem();

  auto problem = problem_builder->ReturnProblem();
  hephaestus::InputParameters exec_params;
  exec_params.SetParam("VisualisationSteps", int(1));
  exec_params.SetParam("UseGLVis", true);
  exec_params.SetParam("Problem", static_cast<hephaestus::SteadyStateProblem *>(problem.get()));

  auto executioner = std::make_unique<hephaestus::SteadyExecutioner>(exec_params);

  hephaestus::logger.info("Created executioner");
  executioner->Execute();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double t;
  double force;
  int block_id;
  for (std::size_t i = 0; i < magnetic_force_monitor->_times.Size(); ++i)
  {
    if (rank == 0)
    {
      t = magnetic_force_monitor->_times[i];
      block_id = magnetic_force_monitor->_attr_out[i];
      force = magnetic_force_monitor->_forces[i];
      hephaestus::logger.info("t = {} s, block {}, F = {} N", t, block_id, force);
    }
  }

  MPI_Finalize();
}