#include "hephaestus.hpp"
#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

class TestEBFormRod
{
protected:
  static double PotentialHigh(const mfem::Vector & x, double t)
  {
    double wj(2.0 * M_PI / 60.0);
    return 2 * cos(wj * t);
  }
  static double PotentialGround(const mfem::Vector & x, double t) { return 0.0; }
  static void EdotBc(const mfem::Vector & x, mfem::Vector & E) { E = 0.0; }
  static void JSrc(const mfem::Vector & x, double t, mfem::Vector & j)
  {
    double wj(2.0 * M_PI / 60.0);
    j[0] = 0.0;
    j[1] = 0.0;
    j[2] = 2 * sin(wj * t);
  }

  hephaestus::InputParameters TestParams()
  {
    double sigma = 2.0 * M_PI * 10;

    double sigma_air;

    sigma_air = 1.0e-6 * sigma;

    hephaestus::Subdomain wire("wire", 1);
    wire.scalar_coefficients.Register(
        "electrical_conductivity", new mfem::ConstantCoefficient(sigma), true);

    hephaestus::Subdomain air("air", 2);
    air.scalar_coefficients.Register(
        "electrical_conductivity", new mfem::ConstantCoefficient(sigma_air), true);

    hephaestus::Coefficients coefficients(std::vector<hephaestus::Subdomain>({wire, air}));

    // coefficients.scalars.Register(
    //     "electrical_conductivity",
    //     new mfem::PWCoefficient(coefficients.getGlobalScalarProperty(
    //         std::string("electrical_conductivity"))),
    //     true);

    hephaestus::BCMap bc_map;
    auto * edot_vec_coef = new mfem::VectorFunctionCoefficient(3, EdotBc);
    bc_map.Register("tangential_dEdt",
                    new hephaestus::VectorDirichletBC(
                        std::string("electric_field"), mfem::Array<int>({1, 2, 3}), edot_vec_coef),
                    true);
    coefficients.scalars.Register(
        "magnetic_permeability", new mfem::ConstantCoefficient(1.0), true);
    coefficients.vectors.Register("surface_tangential_dEdt", edot_vec_coef, true);

    mfem::Array<int> high_terminal(1);
    high_terminal[0] = 1;
    auto * potential_src = new mfem::FunctionCoefficient(PotentialHigh);
    bc_map.Register("high_potential",
                    new hephaestus::ScalarDirichletBC(
                        std::string("electric_potential"), high_terminal, potential_src),
                    true);
    coefficients.scalars.Register("source_potential", potential_src, true);

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 2;
    bc_map.Register(
        "ground_potential",
        new hephaestus::ScalarDirichletBC(std::string("electric_potential"),
                                          ground_terminal,
                                          new mfem::FunctionCoefficient(PotentialGround)),
        true);

    mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./cylinder-hex-q2.gen")).c_str(), 1, 1);

    hephaestus::Outputs outputs;
    outputs.Register("VisItDataCollection", new mfem::VisItDataCollection("EBFormVisIt"), true);
    outputs.Register(
        "ParaViewDataCollection", new mfem::ParaViewDataCollection("EBFormParaView"), true);

    hephaestus::GridFunctions gridfunctions;
    hephaestus::AuxSolvers preprocessors;
    hephaestus::AuxSolvers postprocessors;
    hephaestus::Sources sources;
    hephaestus::InputParameters scalar_potential_source_params;
    scalar_potential_source_params.SetParam("GradPotentialName", std::string("source"));
    scalar_potential_source_params.SetParam("PotentialName", std::string("electric_potential"));
    scalar_potential_source_params.SetParam("HCurlFESpaceName", std::string("HCurl"));
    scalar_potential_source_params.SetParam("H1FESpaceName", std::string("H1"));
    scalar_potential_source_params.SetParam("ConductivityCoefName",
                                            std::string("electrical_conductivity"));
    hephaestus::InputParameters current_solver_options;
    current_solver_options.SetParam("Tolerance", float(1.0e-9));
    current_solver_options.SetParam("MaxIter", (unsigned int)1000);
    current_solver_options.SetParam("PrintLevel", -1);
    scalar_potential_source_params.SetParam("SolverOptions", current_solver_options);
    sources.Register(
        "source", new hephaestus::ScalarPotentialSource(scalar_potential_source_params), true);

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-9));
    solver_options.SetParam("MaxIter", (unsigned int)1000);
    solver_options.SetParam("PrintLevel", 0);

    hephaestus::InputParameters params;
    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("GridFunctions", gridfunctions);
    params.SetParam("PreProcessors", preprocessors);
    params.SetParam("PostProcessors", postprocessors);
    params.SetParam("Sources", sources);
    params.SetParam("Outputs", outputs);
    params.SetParam("SolverOptions", solver_options);

    return params;
  }
};

TEST_CASE_METHOD(TestEBFormRod, "TestEBFormRod", "[CheckRun]")
{
  hephaestus::InputParameters params(TestParams());

  auto problem_builder = std::make_unique<hephaestus::EBDualFormulation>("magnetic_reluctivity",
                                                                         "magnetic_permeability",
                                                                         "electrical_conductivity",
                                                                         "electric_field",
                                                                         "magnetic_flux_density");

  auto bc_map(params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
  auto coefficients(params.GetParam<hephaestus::Coefficients>("Coefficients"));
  auto preprocessors(params.GetParam<hephaestus::AuxSolvers>("PreProcessors"));
  auto postprocessors(params.GetParam<hephaestus::AuxSolvers>("PostProcessors"));
  auto sources(params.GetParam<hephaestus::Sources>("Sources"));
  auto outputs(params.GetParam<hephaestus::Outputs>("Outputs"));
  auto solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
      "SolverOptions", hephaestus::InputParameters()));

  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(params.GetParam<mfem::ParMesh>("Mesh"));
  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P2"));
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P2"));
  problem_builder->AddGridFunction(std::string("analytic_vector_potential"), std::string("HCurl"));
  problem_builder->SetBoundaryConditions(bc_map);
  problem_builder->SetAuxSolvers(preprocessors);
  problem_builder->SetCoefficients(coefficients);
  problem_builder->SetPostprocessors(postprocessors);
  problem_builder->SetSources(sources);
  problem_builder->SetOutputs(outputs);
  problem_builder->SetSolverOptions(solver_options);

  hephaestus::ProblemBuildSequencer sequencer(problem_builder.get());
  sequencer.ConstructEquationSystemProblem();
  std::unique_ptr<hephaestus::TimeDomainProblem> problem = problem_builder->ReturnProblem();
  hephaestus::InputParameters exec_params;
  exec_params.SetParam("TimeStep", float(0.5));
  exec_params.SetParam("StartTime", float(0.00));
  exec_params.SetParam("EndTime", float(2.5));
  exec_params.SetParam("VisualisationSteps", int(1));
  exec_params.SetParam("Problem", problem.get());

  auto executioner = std::make_unique<hephaestus::TransientExecutioner>(exec_params);

  executioner->Execute();
}
