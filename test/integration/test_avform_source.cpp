// Based on an H form MMS test provided by Joseph Dean
#include "hephaestus.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestAVFormSource : public testing::Test {
protected:
  const int var_order{2};
  static double estimate_convergence_rate(HYPRE_BigInt n_i, HYPRE_BigInt n_imo,
                                          double error_i, double error_imo,
                                          int dim) {
    return std::log(error_i / error_imo) /
           std::log(std::pow(n_imo / static_cast<double>(n_i), 1.0 / dim));
  }

  static double potential_ground(const mfem::Vector &x, double t) {
    return -x(0);
  }

  static void adot_bc(const mfem::Vector &x, double t, mfem::Vector &H) {
    H(0) = 1 + sin(x(1) * M_PI) * sin(x(2) * M_PI);
    H(1) = 0;
    H(2) = 0;
  }

  static void A_exact_expr(const mfem::Vector &x, double t,
                           mfem::Vector &A_exact) {
    A_exact(0) = (1 + sin(x(1) * M_PI) * sin(x(2) * M_PI)) * t;
    A_exact(1) = 0;
    A_exact(2) = 0;
  }
  static double mu_expr(const mfem::Vector &x) {
    double variation_scale = 0.0;
    double mu =
        1.0 / (1.0 + variation_scale * cos(M_PI * x(0)) * cos(M_PI * x(1)));
    return mu;
  }

  static void source_field(const mfem::Vector &x, double t, mfem::Vector &f) {
    f(0) = sin(M_PI * x(1)) * sin(M_PI * x(2)) * (t * 2 * M_PI * M_PI + 1);
    f(1) = 0.0;
    f(2) = 0.0;
  }

  hephaestus::InputParameters test_params() {
    hephaestus::Subdomain wire("wire", 1);
    wire.scalar_coefficients.Register("electrical_conductivity",
                                      new mfem::ConstantCoefficient(1.0), true);
    hephaestus::Subdomain air("air", 2);
    air.scalar_coefficients.Register("electrical_conductivity",
                                     new mfem::ConstantCoefficient(1.0), true);

    hephaestus::Coefficients coefficients(
        std::vector<hephaestus::Subdomain>({wire, air}));

    coefficients.scalars.Register("magnetic_permeability",
                                  new mfem::FunctionCoefficient(mu_expr), true);

    hephaestus::BCMap bc_map;
    mfem::VectorFunctionCoefficient *adotVecCoef =
        new mfem::VectorFunctionCoefficient(3, adot_bc);
    bc_map.Register("tangential_dAdt",
                    new hephaestus::VectorDirichletBC(
                        std::string("dmagnetic_vector_potential_dt"),
                        mfem::Array<int>({1, 2, 3}), adotVecCoef),
                    true);
    coefficients.vectors.Register("surface_tangential_dAdt", adotVecCoef, true);
    coefficients.scalars.Register("electrical_conductivity",
                                  new mfem::ConstantCoefficient(1.0), true);

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 1;
    mfem::FunctionCoefficient *ground_coeff =
        new mfem::FunctionCoefficient(potential_ground);
    bc_map.Register("ground_potential",
                    new hephaestus::ScalarDirichletBC(
                        std::string("electric_potential"),
                        mfem::Array<int>({1, 2, 3}), ground_coeff),
                    true);
    coefficients.scalars.Register("ground_potential", ground_coeff, true);

    mfem::VectorFunctionCoefficient *A_exact =
        new mfem::VectorFunctionCoefficient(3, A_exact_expr);
    coefficients.vectors.Register("a_exact_coeff", A_exact, true);

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./beam-tet.mesh")).c_str(), 1, 1);

    hephaestus::Outputs outputs;
    outputs.Register("VisItDataCollection",
                     new mfem::VisItDataCollection("AVFormVisIt"), true);

    hephaestus::InputParameters l2errpostprocparams;
    l2errpostprocparams.SetParam("VariableName",
                                 std::string("magnetic_vector_potential"));
    l2errpostprocparams.SetParam("VectorCoefficientName",
                                 std::string("a_exact_coeff"));
    hephaestus::AuxSolvers postprocessors;
    postprocessors.Register(
        "L2ErrorPostprocessor",
        new hephaestus::L2ErrorVectorPostprocessor(l2errpostprocparams), true);

    hephaestus::VectorCoefficientAux *vec_coef_aux =
        new hephaestus::VectorCoefficientAux("analytic_vector_potential",
                                             "a_exact_coeff");
    vec_coef_aux->SetPriority(-1);
    postprocessors.Register("VectorCoefficientAux", vec_coef_aux, true);

    hephaestus::AuxSolvers preprocessors;

    hephaestus::Sources sources;
    mfem::VectorFunctionCoefficient *JSrcCoef =
        new mfem::VectorFunctionCoefficient(3, source_field);
    coefficients.vectors.Register("source", JSrcCoef, true);
    hephaestus::InputParameters div_free_source_params;
    div_free_source_params.SetParam("SourceName", std::string("source"));
    div_free_source_params.SetParam("HCurlFESpaceName",
                                    std::string("_HCurlFESpace"));
    div_free_source_params.SetParam("H1FESpaceName", std::string("_H1FESpace"));
    div_free_source_params.SetParam("HelmholtzProjection", false);
    sources.Register(
        "source", new hephaestus::DivFreeSource(div_free_source_params), true);

    hephaestus::InputParameters params;
    params.SetParam("TimeStep", float(0.05));
    params.SetParam("StartTime", float(0.00));
    params.SetParam("EndTime", float(0.05));
    params.SetParam("VisualisationSteps", int(1));
    params.SetParam("UseGLVis", true);

    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("PreProcessors", preprocessors);
    params.SetParam("PostProcessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);

    return params;
  }
};

TEST_F(TestAVFormSource, CheckRun) {
  hephaestus::InputParameters params(test_params());
  mfem::ParMesh unrefined_pmesh(params.GetParam<mfem::ParMesh>("Mesh"));

  int num_conv_refinements = 3;
  for (int par_ref_levels = 0; par_ref_levels < num_conv_refinements;
       ++par_ref_levels) {

    std::shared_ptr<mfem::ParMesh> pmesh =
        std::make_shared<mfem::ParMesh>(unrefined_pmesh);

    for (int l = 0; l < par_ref_levels; l++) {
      pmesh->UniformRefinement();
    }
    hephaestus::TimeDomainProblemBuilder *problem_builder =
        new hephaestus::AVFormulation(
            "magnetic_reluctivity", "magnetic_permeability",
            "electrical_conductivity", "magnetic_vector_potential",
            "electric_potential");
    hephaestus::BCMap bc_map(
        params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
    hephaestus::Coefficients coefficients(
        params.GetParam<hephaestus::Coefficients>("Coefficients"));
    //   hephaestus::FESpaces fespaces(
    //       params.GetParam<hephaestus::FESpaces>("FESpaces"));
    //   hephaestus::GridFunctions gridfunctions(
    //       params.GetParam<hephaestus::GridFunctions>("GridFunctions"));
    hephaestus::AuxSolvers preprocessors(
        params.GetParam<hephaestus::AuxSolvers>("PreProcessors"));
    hephaestus::AuxSolvers postprocessors(
        params.GetParam<hephaestus::AuxSolvers>("PostProcessors"));
    hephaestus::Sources sources(
        params.GetParam<hephaestus::Sources>("Sources"));
    hephaestus::Outputs outputs(
        params.GetParam<hephaestus::Outputs>("Outputs"));
    hephaestus::InputParameters solver_options(
        params.GetOptionalParam<hephaestus::InputParameters>(
            "SolverOptions", hephaestus::InputParameters()));

    problem_builder->SetMesh(pmesh);
    problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P2"));
    problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P2"));
    problem_builder->AddGridFunction(std::string("analytic_vector_potential"),
                                     std::string("HCurl"));
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
    exec_params.SetParam("TimeStep", float(0.05));
    exec_params.SetParam("StartTime", float(0.00));
    exec_params.SetParam("EndTime", float(0.05));
    exec_params.SetParam("VisualisationSteps", int(1));
    exec_params.SetParam("UseGLVis", true);
    exec_params.SetParam("Problem", problem.get());
    hephaestus::TransientExecutioner *executioner =
        new hephaestus::TransientExecutioner(exec_params);

    executioner->Init();
    executioner->Execute();
  }

  hephaestus::L2ErrorVectorPostprocessor l2errpostprocessor =
      *(dynamic_cast<hephaestus::L2ErrorVectorPostprocessor *>(
          params.GetParam<hephaestus::AuxSolvers>("PostProcessors")
              .Get("L2ErrorPostprocessor")));

  double r;
  for (std::size_t i = 1; i < l2errpostprocessor.ndofs.Size(); ++i) {
    r = estimate_convergence_rate(
        l2errpostprocessor.ndofs[i], l2errpostprocessor.ndofs[i - 1],
        l2errpostprocessor.l2_errs[i], l2errpostprocessor.l2_errs[i - 1], 3);
    std::cout << r << std::endl;
    ASSERT_TRUE(r > var_order - 0.15);
    ASSERT_TRUE(r < var_order + 1.0);
  }
}
