#include "hephaestus.hpp"
#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

class TestComplexAFormRod
{
protected:
  static double potential_high(const mfem::Vector & x, double t)
  {
    // double wj_(2.0 * M_PI / 60.0);
    // return 2 * cos(wj_ * t);
    return 2.0;
  }

  static double potential_ground(const mfem::Vector & x, double t) { return 0.0; }

  static void a_bc_r(const mfem::Vector & x, mfem::Vector & A)
  {
    A.SetSize(3);
    A = 0.0;
  }

  static void a_bc_i(const mfem::Vector & x, mfem::Vector & A)
  {
    A.SetSize(3);
    A = 0.0;
  }

  hephaestus::InputParameters test_params()
  {
    double sigma = 2.0 * M_PI * 10;

    double sigmaAir;

    sigmaAir = 1.0e-6 * sigma;

    hephaestus::Subdomain wire("wire", 1);
    wire.scalar_coefficients.Register(
        "electrical_conductivity", new mfem::ConstantCoefficient(sigma), true);

    hephaestus::Subdomain air("air", 2);
    air.scalar_coefficients.Register(
        "electrical_conductivity", new mfem::ConstantCoefficient(sigmaAir), true);

    hephaestus::Coefficients coefficients(std::vector<hephaestus::Subdomain>({wire, air}));

    coefficients.scalars.Register("frequency", new mfem::ConstantCoefficient(1.0 / 60.0), true);
    coefficients.scalars.Register(
        "dielectric_permittivity", new mfem::ConstantCoefficient(0.0), true);
    coefficients.scalars.Register(
        "magnetic_permeability", new mfem::ConstantCoefficient(1.0), true);

    hephaestus::BCMap bc_map;

    bc_map.Register(
        "tangential_A",
        new hephaestus::VectorDirichletBC(std::string("magnetic_vector_potential"),
                                          mfem::Array<int>({1, 2, 3}),
                                          new mfem::VectorFunctionCoefficient(3, a_bc_r),
                                          new mfem::VectorFunctionCoefficient(3, a_bc_i)),
        true);

    mfem::Array<int> high_terminal(1);
    high_terminal[0] = 1;
    mfem::FunctionCoefficient * potential_src = new mfem::FunctionCoefficient(potential_high);
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
                                          new mfem::FunctionCoefficient(potential_ground)),
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
    scalar_potential_source_params.SetParam("SourceName", std::string("source"));
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
    params.SetParam("UseGLVis", true);
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

TEST_CASE_METHOD(TestComplexAFormRod, "TestComplexAFormRod", "[CheckRun]")
{
  hephaestus::InputParameters params(test_params());
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(params.GetParam<mfem::ParMesh>("Mesh"));
  hephaestus::ComplexAFormulation * problem_builder =
      new hephaestus::ComplexAFormulation("magnetic_reluctivity",
                                          "electrical_conductivity",
                                          "dielectric_permittivity",
                                          "frequency",
                                          "magnetic_vector_potential",
                                          "magnetic_vector_potential_real",
                                          "magnetic_vector_potential_imag");
  hephaestus::BCMap bc_map(params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
  hephaestus::Coefficients coefficients(params.GetParam<hephaestus::Coefficients>("Coefficients"));
  hephaestus::AuxSolvers preprocessors(params.GetParam<hephaestus::AuxSolvers>("PreProcessors"));
  hephaestus::AuxSolvers postprocessors(params.GetParam<hephaestus::AuxSolvers>("PostProcessors"));
  hephaestus::Sources sources(params.GetParam<hephaestus::Sources>("Sources"));
  hephaestus::Outputs outputs(params.GetParam<hephaestus::Outputs>("Outputs"));
  hephaestus::InputParameters solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
      "SolverOptions", hephaestus::InputParameters()));

  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P1"));
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P1"));
  problem_builder->AddGridFunction("magnetic_vector_potential_real", "HCurl");
  problem_builder->AddGridFunction("magnetic_vector_potential_imag", "HCurl");
  problem_builder->SetBoundaryConditions(bc_map);
  problem_builder->SetAuxSolvers(preprocessors);
  problem_builder->SetCoefficients(coefficients);
  problem_builder->SetPostprocessors(postprocessors);

  problem_builder->AddFESpace("HDiv", "RT_3D_P0");
  problem_builder->AddFESpace("Scalar_L2", "L2Int_3D_P0");

  problem_builder->AddGridFunction("electric_field_real", "HCurl");
  problem_builder->AddGridFunction("electric_field_imag", "HCurl");
  problem_builder->registerElectricFieldAux("electric_field_real", "electric_field_imag");

  problem_builder->AddGridFunction("magnetic_flux_density_real", "HDiv");
  problem_builder->AddGridFunction("magnetic_flux_density_imag", "HDiv");
  problem_builder->registerMagneticFluxDensityAux("magnetic_flux_density_real",
                                                  "magnetic_flux_density_imag");

  problem_builder->AddGridFunction("current_density_real", "HDiv");
  problem_builder->AddGridFunction("current_density_imag", "HDiv");
  problem_builder->registerCurrentDensityAux("current_density_real", "current_density_imag");

  problem_builder->AddGridFunction("joule_heating_density", "Scalar_L2");
  problem_builder->registerJouleHeatingDensityAux("joule_heating_density",
                                                  "electric_field_real",
                                                  "electric_field_imag",
                                                  "electrical_conductivity");

  problem_builder->SetSources(sources);
  problem_builder->SetOutputs(outputs);
  problem_builder->SetSolverOptions(solver_options);

  hephaestus::ProblemBuildSequencer sequencer(problem_builder);
  sequencer.ConstructOperatorProblem();
  std::unique_ptr<hephaestus::SteadyStateProblem> problem = problem_builder->ReturnProblem();

  hephaestus::InputParameters exec_params;
  exec_params.SetParam("Problem", problem.get());
  hephaestus::SteadyExecutioner * executioner = new hephaestus::SteadyExecutioner(exec_params);

  executioner->Execute();
}
