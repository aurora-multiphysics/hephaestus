// Based on an H form MMS test provided by Joseph Dean
#include "hephaestus.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestAFormSource : public testing::Test {
protected:
  static double estimate_convergence_rate(HYPRE_BigInt n_i, HYPRE_BigInt n_imo,
                                          double error_i, double error_imo,
                                          int dim) {
    return std::log(error_i / error_imo) /
           std::log(std::pow(n_imo / static_cast<double>(n_i), 1.0 / dim));
  }

  static double potential_ground(const mfem::Vector &x, double t) {
    return 0.0;
  }

  static void adot_bc(const mfem::Vector &x, double t, mfem::Vector &A) {
    A(0) = sin(x(1) * M_PI) * sin(x(2) * M_PI);
    A(1) = 0;
    A(2) = 0;
  }

  static void A_exact_expr(const mfem::Vector &x, double t,
                           mfem::Vector &A_exact) {
    A_exact(0) = sin(x(1) * M_PI) * sin(x(2) * M_PI) * t;
    A_exact(1) = 0;
    A_exact(2) = 0;
  }
  static double mu_expr(const mfem::Vector &x) {
    double variation_scale = 0.0;
    double mu =
        1.0 / (1.0 + variation_scale * cos(M_PI * x(0)) * cos(M_PI * x(1)));
    return mu;
  }
  // Source field
  static void source_field(const mfem::Vector &x, double t, mfem::Vector &f) {
    double variation_scale = 0.0;
    f(0) = t * M_PI * M_PI * sin(M_PI * x(1)) * sin(M_PI * x(2)) *
               (3 * variation_scale * cos(M_PI * x(0)) * cos(M_PI * x(1)) + 2) +
           sin(M_PI * x(1)) * sin(M_PI * x(2));
    f(1) = -variation_scale * M_PI * M_PI * t * sin(M_PI * x(0)) *
           cos(M_PI * x(1)) * cos(M_PI * x(1)) * sin(M_PI * x(2));
    f(2) = -0.5 * variation_scale * M_PI * M_PI * t * sin(M_PI * x(0)) *
           sin(2 * M_PI * x(1)) * cos(M_PI * x(2));
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
                    new hephaestus::VectorFunctionDirichletBC(
                        std::string("dmagnetic_vector_potential_dt"),
                        mfem::Array<int>({1, 2, 3}), adotVecCoef),
                    true);
    coefficients.vectors.Register("surface_tangential_dAdt", adotVecCoef, true);
    coefficients.scalars.Register("electrical_conductivity",
                                  new mfem::ConstantCoefficient(1.0), true);

    mfem::VectorFunctionCoefficient *A_exact =
        new mfem::VectorFunctionCoefficient(3, A_exact_expr);
    coefficients.vectors.Register("a_exact_coeff", A_exact, true);

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./beam-tet.mesh")).c_str(), 1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("AFormVisIt");
    hephaestus::Outputs outputs(data_collections);

    hephaestus::InputParameters l2errpostprocparams;
    l2errpostprocparams.SetParam("VariableName",
                                 std::string("magnetic_vector_potential"));
    l2errpostprocparams.SetParam("VectorCoefficientName",
                                 std::string("a_exact_coeff"));
    hephaestus::AuxSolvers postprocessors;
    postprocessors.Register(
        "L2ErrorPostprocessor",
        new hephaestus::L2ErrorVectorPostprocessor(l2errpostprocparams), true);

    hephaestus::InputParameters vectorcoeffauxparams;
    vectorcoeffauxparams.SetParam("VariableName",
                                  std::string("analytic_vector_potential"));
    vectorcoeffauxparams.SetParam("VectorCoefficientName",
                                  std::string("a_exact_coeff"));

    hephaestus::VectorCoefficientAuxSolver *vec_coef_aux =
        new hephaestus::VectorCoefficientAuxSolver(vectorcoeffauxparams);
    vec_coef_aux->SetPriority(-1);
    postprocessors.Register("VectorCoefficientAuxSolver", vec_coef_aux, true);

    hephaestus::Sources sources;
    mfem::VectorFunctionCoefficient *JSrcCoef =
        new mfem::VectorFunctionCoefficient(3, source_field);
    coefficients.vectors.Register("source", JSrcCoef, true);
    hephaestus::InputParameters div_free_source_params;
    div_free_source_params.SetParam("SourceName", std::string("source"));
    div_free_source_params.SetParam("HCurlFESpaceName",
                                    std::string("_HCurlFESpace"));
    div_free_source_params.SetParam("H1FESpaceName", std::string("H1"));
    hephaestus::InputParameters current_solver_options;
    current_solver_options.SetParam("Tolerance", float(1.0e-12));
    current_solver_options.SetParam("MaxIter", (unsigned int)200);
    current_solver_options.SetParam("PrintLevel", 0);
    div_free_source_params.SetParam("SolverOptions", current_solver_options);
    div_free_source_params.SetParam("HelmholtzProjection", false);
    sources.Register(
        "source", new hephaestus::DivFreeSource(div_free_source_params), true);

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-16));
    solver_options.SetParam("MaxIter", (unsigned int)1000);
    solver_options.SetParam("PrintLevel", 0);

    hephaestus::InputParameters params;
    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("PostProcessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);
    params.SetParam("SolverOptions", solver_options);

    return params;
  }
};

TEST_F(TestAFormSource, CheckRun) {
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
        new hephaestus::AFormulation();
    hephaestus::BCMap bc_map(
        params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
    hephaestus::Coefficients coefficients(
        params.GetParam<hephaestus::Coefficients>("Coefficients"));
    //   hephaestus::FESpaces fespaces(
    //       params.GetParam<hephaestus::FESpaces>("FESpaces"));
    //   hephaestus::GridFunctions gridfunctions(
    //       params.GetParam<hephaestus::GridFunctions>("GridFunctions"));
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

    delete executioner;
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
    ASSERT_TRUE(r > 2 - 0.15);
    ASSERT_TRUE(r < 2 + 1.0);
  }
}
