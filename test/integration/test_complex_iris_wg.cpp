#include "hephaestus.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

extern const char *DATA_DIR;

class TestComplexIrisWaveguide {
protected:
  static void e_bc_r(const mfem::Vector &x, mfem::Vector &E) {
    E.SetSize(3);
    E = 0.0;
  }

  static void e_bc_i(const mfem::Vector &x, mfem::Vector &E) {
    E.SetSize(3);
    E = 0.0;
  }

  inline static const double epsilon0_ = 8.8541878176e-12; // F/m;
  inline static const double mu0_ = 4.0e-7 * M_PI;         // H/m;
  inline static const double freq_ = 9.3e9;                // 10/2pi
  double port_length_vector[3] = {0.0, 22.86e-3, 0.0};
  double port_width_vector[3] = {0.0, 0.0, 10.16e-3};

  hephaestus::InputParameters test_params() {
    hephaestus::Subdomain air("air", 1);

    air.scalar_coefficients.Register("real_electrical_conductivity",
                                     new mfem::ConstantCoefficient(0.0), true);
    air.scalar_coefficients.Register("imag_electrical_conductivity",
                                     new mfem::ConstantCoefficient(0.0), true);
    air.scalar_coefficients.Register("real_rel_permittivity",
                                     new mfem::ConstantCoefficient(1.0), true);
    air.scalar_coefficients.Register("imag_rel_permittivity",
                                     new mfem::ConstantCoefficient(0.0), true);
    air.scalar_coefficients.Register("real_rel_permeability",
                                     new mfem::ConstantCoefficient(1.0), true);
    air.scalar_coefficients.Register("imag_rel_permeability",
                                     new mfem::ConstantCoefficient(0.0), true);

    hephaestus::Coefficients coefficients(
        std::vector<hephaestus::Subdomain>({air}));

    coefficients.scalars.Register("frequency",
                                  new mfem::ConstantCoefficient(freq_), true);
    coefficients.scalars.Register("magnetic_permeability",
                                  new mfem::ConstantCoefficient(mu0_), true);
    coefficients.scalars.Register("dielectric_permittivity",
                                  new mfem::ConstantCoefficient(epsilon0_),
                                  true);
    coefficients.scalars.Register("electrical_conductivity",
                                  new mfem::ConstantCoefficient(0.0), true);

    hephaestus::BCMap bc_map;
    mfem::Array<int> dirichlet_attr(1);
    dirichlet_attr[0] = 1;
    bc_map.Register("tangential_E",
                    new hephaestus::VectorDirichletBC(
                        std::string("electric_field"), dirichlet_attr,
                        new mfem::VectorFunctionCoefficient(3, e_bc_r),
                        new mfem::VectorFunctionCoefficient(3, e_bc_i)),
                    true);

    mfem::Array<int> wgi_in_attr(1);
    wgi_in_attr[0] = 2;
    bc_map.Register("WaveguidePortIn",
                    new hephaestus::RWTE10PortRBC(
                        std::string("electric_field"), wgi_in_attr, freq_,
                        port_length_vector, port_width_vector, true),
                    true);

    mfem::Array<int> wgi_out_attr(1);
    wgi_out_attr[0] = 3;
    bc_map.Register("WaveguidePortOut",
                    new hephaestus::RWTE10PortRBC(
                        std::string("electric_field"), wgi_out_attr, freq_,
                        port_length_vector, port_width_vector, false),
                    true);

    mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./irises.g")).c_str(),
                    1, 1);

    hephaestus::Outputs outputs;
    outputs.Register("VisItDataCollection",
                     new mfem::VisItDataCollection("Hertz-AMR-Parallel-VisIt"),
                     true);

    hephaestus::FESpaces fespaces;
    hephaestus::GridFunctions gridfunctions;
    hephaestus::AuxSolvers postprocessors;
    hephaestus::AuxSolvers preprocessors;
    hephaestus::Sources sources;

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-16));
    solver_options.SetParam("MaxIter", (unsigned int)1000);
    solver_options.SetParam("PrintLevel", 0);

    hephaestus::InputParameters params;
    params.SetParam("UseGLVis", true);

    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("FESpaces", fespaces);
    params.SetParam("GridFunctions", gridfunctions);
    params.SetParam("PreProcessors", preprocessors);
    params.SetParam("PostProcessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);
    params.SetParam("SolverOptions", solver_options);

    return params;
  }
};

TEST_CASE_METHOD(TestComplexIrisWaveguide, "TestComplexIrisWaveguide",
                 "[CheckRun]") {
  hephaestus::InputParameters params(test_params());
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(params.GetParam<mfem::ParMesh>("Mesh"));

  auto problem_builder = new hephaestus::ComplexEFormulation(
      "magnetic_reluctivity", "electrical_conductivity",
      "dielectric_permittivity", "frequency", "electric_field",
      "electric_field_real", "electric_field_imag");

  hephaestus::BCMap bc_map(
      params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
  hephaestus::Coefficients coefficients(
      params.GetParam<hephaestus::Coefficients>("Coefficients"));
  hephaestus::AuxSolvers preprocessors(
      params.GetParam<hephaestus::AuxSolvers>("PreProcessors"));
  hephaestus::AuxSolvers postprocessors(
      params.GetParam<hephaestus::AuxSolvers>("PostProcessors"));
  hephaestus::Sources sources(params.GetParam<hephaestus::Sources>("Sources"));
  hephaestus::Outputs outputs(params.GetParam<hephaestus::Outputs>("Outputs"));
  hephaestus::InputParameters solver_options(
      params.GetOptionalParam<hephaestus::InputParameters>(
          "SolverOptions", hephaestus::InputParameters()));

  problem_builder->SetMesh(pmesh);
  problem_builder->SetBoundaryConditions(bc_map);
  problem_builder->SetAuxSolvers(preprocessors);

  problem_builder->AddFESpace("HDiv", "RT_3D_P0");
  problem_builder->AddFESpace("Scalar_L2", "L2Int_3D_P0");

  problem_builder->AddGridFunction("magnetic_flux_density_real", "HDiv");
  problem_builder->AddGridFunction("magnetic_flux_density_imag", "HDiv");
  problem_builder->registerMagneticFluxDensityAux("magnetic_flux_density_real",
                                                  "magnetic_flux_density_imag");

  problem_builder->AddGridFunction("current_density_real", "HDiv");
  problem_builder->AddGridFunction("current_density_imag", "HDiv");
  problem_builder->registerCurrentDensityAux("current_density_real",
                                             "current_density_imag");

  problem_builder->AddGridFunction("joule_heating_density", "Scalar_L2");
  problem_builder->registerJouleHeatingDensityAux(
      "joule_heating_density", "electric_field_real", "electric_field_imag",
      "electrical_conductivity");

  problem_builder->SetCoefficients(coefficients);
  problem_builder->SetPostprocessors(postprocessors);
  problem_builder->SetSources(sources);
  problem_builder->SetOutputs(outputs);
  problem_builder->SetSolverOptions(solver_options);

  hephaestus::ProblemBuildSequencer sequencer(problem_builder);
  sequencer.ConstructOperatorProblem();
  std::unique_ptr<hephaestus::SteadyStateProblem> problem =
      problem_builder->ReturnProblem();

  hephaestus::InputParameters exec_params;
  exec_params.SetParam("Problem", problem.get());

  auto executioner =
      std::make_unique<hephaestus::SteadyExecutioner>(exec_params);

  executioner->Execute();

  mfem::Vector zeroVec(3);
  zeroVec = 0.0;
  mfem::VectorConstantCoefficient zeroCoef(zeroVec);

  double norm_r =
      executioner->problem->gridfunctions.Get("electric_field_real")
          ->ComputeMaxError(zeroCoef);
  double norm_i =
      executioner->problem->gridfunctions.Get("electric_field_imag")
          ->ComputeMaxError(zeroCoef);
  REQUIRE_THAT(norm_r, Catch::Matchers::WithinAbs(4896.771, 0.001));
  REQUIRE_THAT(norm_i, Catch::Matchers::WithinAbs(5357.650, 0.001));
}
