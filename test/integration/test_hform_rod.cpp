#include "hephaestus.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestHFormRod : public testing::Test {
protected:
  static double potential_high(const mfem::Vector &x, double t) {
    double wj_(2.0 * M_PI / 60.0);
    return cos(wj_ * t);
  }
  static double potential_ground(const mfem::Vector &x, double t) {
    return 0.0;
  }
  static void hdot_bc(const mfem::Vector &x, mfem::Vector &H) { H = 0.0; }
  static void b_src(const mfem::Vector &x, double t, mfem::Vector &b) {
    double wj_(2.0 * M_PI / 60.0);
    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = 2 * sin(wj_ * t);
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

    // coefficients.scalars.Register(
    //     "electrical_conductivity",
    //     new mfem::PWCoefficient(coefficients.getGlobalScalarProperty(
    //         std::string("electrical_conductivity"))),
    //     true);

    hephaestus::BCMap bc_map;
    mfem::VectorFunctionCoefficient *hdotVecCoef =
        new mfem::VectorFunctionCoefficient(3, hdot_bc);
    bc_map.Register("tangential_dHdt",
                    new hephaestus::VectorDirichletBC(
                        std::string("dmagnetic_field_dt"),
                        mfem::Array<int>({1, 2, 3}), hdotVecCoef),
                    true);
    coefficients.vectors.Register("surface_tangential_dHdt", hdotVecCoef, true);
    coefficients.scalars.Register("magnetic_permeability",
                                  new mfem::ConstantCoefficient(1.0), true);

    mfem::Array<int> high_terminal(1);
    high_terminal[0] = 1;
    bc_map.Register("high_potential",
                    new hephaestus::ScalarDirichletBC(
                        std::string("magnetic_potential"), high_terminal,
                        new mfem::FunctionCoefficient(potential_high)),
                    true);

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 2;
    bc_map.Register("ground_potential",
                    new hephaestus::ScalarDirichletBC(
                        std::string("magnetic_potential"), ground_terminal,
                        new mfem::FunctionCoefficient(potential_ground)),
                    true);

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./cylinder-hex-q2.gen")).c_str(),
        1, 1);

    hephaestus::Outputs outputs;
    outputs.Register("VisItDataCollection",
                     new mfem::VisItDataCollection("HFormVisIt"), true);
    outputs.Register("ParaViewDataCollection",
                     new mfem::ParaViewDataCollection("HFormParaView"), true);

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-16));
    solver_options.SetParam("MaxIter", (unsigned int)1000);
    solver_options.SetParam("PrintLevel", 0);

    hephaestus::GridFunctions gridfunctions;
    hephaestus::AuxSolvers preprocessors;
    hephaestus::AuxSolvers postprocessors;
    hephaestus::Sources sources;
    hephaestus::InputParameters scalar_potential_source_params;
    scalar_potential_source_params.SetParam("EFieldName",
                                            std::string("source"));
    scalar_potential_source_params.SetParam("PotentialName",
                                            std::string("magnetic_potential"));
    scalar_potential_source_params.SetParam("HCurlFESpaceName",
                                            std::string("_HCurlFESpace"));
    scalar_potential_source_params.SetParam("H1FESpaceName", std::string("H1"));
    scalar_potential_source_params.SetParam(
        "ConductivityCoefName", std::string("magnetic_permeability"));
    hephaestus::InputParameters current_solver_options;
    current_solver_options.SetParam("Tolerance", float(1.0e-12));
    current_solver_options.SetParam("MaxIter", (unsigned int)200);
    current_solver_options.SetParam("PrintLevel", 0);
    scalar_potential_source_params.SetParam("SolverOptions",
                                            current_solver_options);
    sources.Register(
        "source",
        new hephaestus::ScalarPotentialSource(scalar_potential_source_params),
        true);

    hephaestus::InputParameters params;
    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("GridFunctions", gridfunctions);
    params.SetParam("PreProcessors", preprocessors);
    params.SetParam("PostProcessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);
    params.SetParam("SolverOptions", solver_options);

    return params;
  }
};

TEST_F(TestHFormRod, CheckRun) {
  hephaestus::InputParameters params(test_params());
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(params.GetParam<mfem::ParMesh>("Mesh"));

  hephaestus::TimeDomainProblemBuilder *problem_builder =
      new hephaestus::HFormulation("electrical_resistivity",
                                   "electrical_conductivity",
                                   "magnetic_permeability", "magnetic_field");
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
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P2"));
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
