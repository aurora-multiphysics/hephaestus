#include "hephaestus.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestAVFormRod : public testing::Test {
protected:
  static double potential_high(const mfem::Vector &x, double t) {
    double wj_(2.0 * M_PI / 60.0);
    return 2 * cos(wj_ * t);
  }
  static double potential_ground(const mfem::Vector &x, double t) {
    return 0.0;
  }

  static void adot_bc(const mfem::Vector &x, double t, mfem::Vector &dAdt) {
    dAdt = 0.0;
  }

  static void source_current(const mfem::Vector &x, double t, mfem::Vector &J) {
    J = 0.0;
  }

  hephaestus::InputParameters test_params() {
    double sigma = 2.0 * M_PI * 10;

    double sigmaAir;

    sigmaAir = 1.0e-6 * sigma;

    hephaestus::Subdomain wire("wire", 1);
    wire.scalar_coefficients.Register(
        "electrical_conductivity", new mfem::ConstantCoefficient(sigma), true);

    hephaestus::Subdomain air("air", 2);
    air.scalar_coefficients.Register("electrical_conductivity",
                                     new mfem::ConstantCoefficient(sigmaAir),
                                     true);

    hephaestus::Coefficients coefficients(
        std::vector<hephaestus::Subdomain>({wire, air}));

    coefficients.scalars.Register("magnetic_permeability",
                                  new mfem::ConstantCoefficient(1.0), true);

    hephaestus::BCMap bc_map;
    mfem::VectorFunctionCoefficient *adotVecCoef =
        new mfem::VectorFunctionCoefficient(3, adot_bc);
    bc_map.Register("tangential_dAdt",
                    new hephaestus::VectorDirichletBC(
                        std::string("dmagnetic_vector_potential_dt"),
                        mfem::Array<int>({1, 2, 3}), adotVecCoef),
                    true);
    coefficients.vectors.Register("surface_tangential_dAdt", adotVecCoef, true);

    mfem::Array<int> high_terminal(1);
    high_terminal[0] = 1;
    bc_map.Register("high_potential",
                    new hephaestus::ScalarDirichletBC(
                        std::string("electric_potential"), high_terminal,
                        new mfem::FunctionCoefficient(potential_high)),
                    true);

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 2;
    bc_map.Register("ground_potential",
                    new hephaestus::ScalarDirichletBC(
                        std::string("electric_potential"), ground_terminal,
                        new mfem::FunctionCoefficient(potential_ground)),
                    true);

    mfem::VectorFunctionCoefficient *JSrcCoef =
        new mfem::VectorFunctionCoefficient(3, source_current);

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./cylinder-hex-q2.gen")).c_str(),
        1, 1);

    hephaestus::Outputs outputs;
    outputs.Register("VisItDataCollection",
                     new mfem::VisItDataCollection("AVFormVisIt"), true);
    outputs.Register("ParaViewDataCollection",
                     new mfem::ParaViewDataCollection("AVFormParaView"), true);

    hephaestus::FESpaces fespaces;
    hephaestus::GridFunctions gridfunctions;
    hephaestus::AuxSolvers postprocessors;
    hephaestus::AuxSolvers preprocessors;
    hephaestus::Sources sources;

    hephaestus::InputParameters params;
    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("FESpaces", fespaces);
    params.SetParam("GridFunctions", gridfunctions);
    params.SetParam("PreProcessors", preprocessors);
    params.SetParam("PostProcessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);

    return params;
  }
};

TEST_F(TestAVFormRod, CheckRun) {
  hephaestus::InputParameters params(test_params());
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(params.GetParam<mfem::ParMesh>("Mesh"));

  hephaestus::TimeDomainProblemBuilder *problem_builder =
      new hephaestus::AVFormulation(
          "magnetic_reluctivity", "magnetic_permeability",
          "electrical_conductivity", "magnetic_vector_potential",
          "electric_potential");
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
  sequencer.ConstructEquationSystemProblem();
  std::unique_ptr<hephaestus::TimeDomainProblem> problem =
      problem_builder->ReturnProblem();

  hephaestus::InputParameters exec_params;
  exec_params.SetParam("TimeStep", float(0.5));
  exec_params.SetParam("StartTime", float(0.00));
  exec_params.SetParam("EndTime", float(2.5));
  exec_params.SetParam("VisualisationSteps", int(1));
  exec_params.SetParam("Problem", problem.get());
  hephaestus::TransientExecutioner *executioner =
      new hephaestus::TransientExecutioner(exec_params);

  executioner->Execute();
}
