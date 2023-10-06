#include "../common/pfem_extras.hpp"
#include "factory.hpp"
#include "inputs.hpp"
#include "steady_executioner.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestComplexERMESMouse : public testing::Test {
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
  inline static const double freq_ = 900e6;                // 10/2pi
  double port_length_vector[3] = {24.76e-2, 0.0, 0.0};
  double port_width_vector[3] = {0.0, 12.38e-2, 0.0};
  //   double port_length_vector[3] = {0.0, 22.86e-3, 0.0};
  //   double port_width_vector[3] = {0.0, 0.0, 10.16e-3};

  hephaestus::InputParameters test_params() {
    hephaestus::Subdomain mouse("mouse", 1);
    hephaestus::Subdomain air("air", 2);

    air.scalar_coefficients.Register("electrical_conductivity",
                                     new mfem::ConstantCoefficient(0.0), true);
    air.scalar_coefficients.Register("dielectric_permittivity",
                                     new mfem::ConstantCoefficient(epsilon0_),
                                     true);
    air.scalar_coefficients.Register("magnetic_permeability",
                                     new mfem::ConstantCoefficient(mu0_), true);

    mouse.scalar_coefficients.Register(
        "electrical_conductivity", new mfem::ConstantCoefficient(0.97), true);
    mouse.scalar_coefficients.Register(
        "dielectric_permittivity",
        new mfem::ConstantCoefficient(43 * epsilon0_), true);
    mouse.scalar_coefficients.Register(
        "magnetic_permeability", new mfem::ConstantCoefficient(mu0_), true);

    hephaestus::Coefficients coefficients(
        std::vector<hephaestus::Subdomain>({air, mouse}));

    coefficients.scalars.Register("frequency",
                                  new mfem::ConstantCoefficient(freq_), true);

    hephaestus::BCMap bc_map;
    mfem::Array<int> dirichlet_attr({2, 3, 4});
    bc_map.Register("tangential_E",
                    new hephaestus::VectorFunctionDirichletBC(
                        std::string("electric_field"), dirichlet_attr,
                        new mfem::VectorFunctionCoefficient(3, e_bc_r),
                        new mfem::VectorFunctionCoefficient(3, e_bc_i)),
                    true);

    mfem::Array<int> wgi_in_attr(1);
    wgi_in_attr[0] = 5;
    bc_map.Register("WaveguidePortIn",
                    new hephaestus::RWTE10PortRBC(
                        std::string("electric_field"), wgi_in_attr, freq_,
                        port_length_vector, port_width_vector, true),
                    true);

    mfem::Array<int> wgi_out_attr(1);
    wgi_out_attr[0] = 6;
    bc_map.Register("WaveguidePortOut",
                    new hephaestus::RWTE10PortRBC(
                        std::string("electric_field"), wgi_out_attr, freq_,
                        port_length_vector, port_width_vector, false),
                    true);

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./ermes_mouse_coarse.g")).c_str(),
        1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("ComplexMaxwellERMESMouse");
    hephaestus::Outputs outputs(data_collections);

    hephaestus::FESpaces fespaces;
    hephaestus::GridFunctions gridfunctions;
    hephaestus::AuxSolvers postprocessors;
    hephaestus::AuxSolvers preprocessors;
    hephaestus::Sources sources;

    hephaestus::FrequencyDomainFormulation *formulation =
        new hephaestus::ComplexEFormulation(
            "magnetic_reluctivity", "electrical_conductivity",
            "dielectric_permittivity", "frequency", "electric_field");

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

TEST_F(TestComplexERMESMouse, CheckRun) {
  hephaestus::InputParameters params(test_params());
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(params.GetParam<mfem::ParMesh>("Mesh"));

  hephaestus::SteadyStateProblemBuilder *problem_builder =
      new hephaestus::ComplexEFormulation(
          "magnetic_reluctivity", "electrical_conductivity",
          "dielectric_permittivity", "frequency", "electric_field");
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
  exec_params.SetParam("UseGLVis", true);
  exec_params.SetParam("Problem", problem.get());
  hephaestus::SteadyExecutioner *executioner =
      new hephaestus::SteadyExecutioner(exec_params);
  executioner->Init();
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
  ASSERT_NEAR(norm_r, 480, 15);
  ASSERT_NEAR(norm_i, 180, 5);
}
