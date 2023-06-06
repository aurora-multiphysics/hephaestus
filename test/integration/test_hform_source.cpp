// Based on an H form MMS test provided by Joseph Dean
#include "hephaestus.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestHFormSource : public testing::Test {
protected:
  const int var_order{2};
  static double estimate_convergence_rate(HYPRE_BigInt n_i, HYPRE_BigInt n_imo,
                                          double error_i, double error_imo,
                                          int dim) {
    return std::log(error_i / error_imo) /
           std::log(std::pow(n_imo / static_cast<double>(n_i), 1.0 / dim));
  }

  static double potential_ground(const mfem::Vector &x, double t) {
    return 0.0;
  }

  static void hdot_bc(const mfem::Vector &x, double t, mfem::Vector &H) {
    H(0) = sin(x(1) * M_PI) * sin(x(2) * M_PI);
    H(1) = 0;
    H(2) = 0;
  }

  static void H_exact_expr(const mfem::Vector &x, double t,
                           mfem::Vector &H_exact) {
    H_exact(0) = sin(x(1) * M_PI) * sin(x(2) * M_PI) * t;
    H_exact(1) = 0;
    H_exact(2) = 0;
  }
  static double sigma_expr(const mfem::Vector &x) {
    double variation_scale = 0.5;
    double sigma =
        1.0 / (1.0 + variation_scale * cos(M_PI * x(0)) * cos(M_PI * x(1)));
    return sigma;
  }
  // Source field
  static void source_field(const mfem::Vector &x, double t, mfem::Vector &f) {
    double variation_scale = 0.5;
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
    wire.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::Subdomain air("air", 2);
    air.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);

    hephaestus::DomainProperties domain_properties(
        std::vector<hephaestus::Subdomain>({wire, air}));

    domain_properties.scalar_property_map["electrical_conductivity"] =
        new mfem::FunctionCoefficient(sigma_expr);

    hephaestus::BCMap bc_map;
    mfem::VectorFunctionCoefficient *hdotVecCoef =
        new mfem::VectorFunctionCoefficient(3, hdot_bc);
    bc_map.Register("tangential_dHdt",
                    new hephaestus::VectorFunctionDirichletBC(
                        std::string("dmagnetic_field_dt"),
                        mfem::Array<int>({1, 2, 3}), hdotVecCoef),
                    true);
    domain_properties.vector_property_map["surface_tangential_dHdt"] =
        hdotVecCoef;
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::ConstantCoefficient(1.0);

    mfem::VectorFunctionCoefficient *dBdtSrcCoef =
        new mfem::VectorFunctionCoefficient(3, source_field);
    domain_properties.vector_property_map["source"] = dBdtSrcCoef;

    // bc_map.Register("ground_potential",
    //                 new hephaestus::FunctionDirichletBC(
    //                     std::string("magnetic_potential"),
    //                     mfem::Array<int>({1, 2, 3}),
    //                     new mfem::FunctionCoefficient(potential_ground)),
    //                 true);
    // bc_map.Register("ground_potential",
    //                 new hephaestus::IntegratedBC(
    //                     std::string("magnetic_potential"),
    //                     mfem::Array<int>({1, 2, 3}),
    //                     new mfem::BoundaryNormalLFIntegrator(*dBdtSrcCoef)),
    //                 true);

    mfem::VectorFunctionCoefficient *H_exact =
        new mfem::VectorFunctionCoefficient(3, H_exact_expr);
    domain_properties.vector_property_map["h_exact_coeff"] = H_exact;

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./beam-tet.mesh")).c_str(), 1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("HFormVisIt");
    hephaestus::Outputs outputs(data_collections);

    hephaestus::InputParameters l2errpostprocparams;
    l2errpostprocparams.SetParam("VariableName", std::string("magnetic_field"));
    l2errpostprocparams.SetParam("VectorCoefficientName",
                                 std::string("h_exact_coeff"));
    hephaestus::AuxSolvers postprocessors;
    postprocessors.Register(
        "L2ErrorPostprocessor",
        new hephaestus::L2ErrorVectorPostprocessor(l2errpostprocparams), true);

    hephaestus::InputParameters vectorcoeffauxparams;
    vectorcoeffauxparams.SetParam("VariableName",
                                  std::string("analytic_magnetic_field"));
    vectorcoeffauxparams.SetParam("VectorCoefficientName",
                                  std::string("h_exact_coeff"));

    hephaestus::AuxSolvers auxsolvers;
    auxsolvers.Register(
        "VectorCoefficientAuxSolver",
        new hephaestus::VectorCoefficientAuxSolver(vectorcoeffauxparams), true);

    hephaestus::Sources sources;
    hephaestus::InputParameters div_free_source_params;
    div_free_source_params.SetParam("SourceName", std::string("source"));
    div_free_source_params.SetParam("PotentialName",
                                    std::string("magnetic_potential"));
    div_free_source_params.SetParam("HCurlFESpaceName",
                                    std::string("_HCurlFESpace"));
    div_free_source_params.SetParam("ConductivityCoefName",
                                    std::string("electrical_conductivity"));
    div_free_source_params.SetParam("H1FESpaceName", std::string("H1"));
    hephaestus::InputParameters current_solver_options;
    current_solver_options.SetParam("Tolerance", float(1.0e-12));
    current_solver_options.SetParam("MaxIter", (unsigned int)200);
    current_solver_options.SetParam("PrintLevel", 0);
    div_free_source_params.SetParam("SolverOptions", current_solver_options);
    sources.Register(
        "source", new hephaestus::DivFreeSource(div_free_source_params), true);

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-16));
    solver_options.SetParam("MaxIter", (unsigned int)1000);
    solver_options.SetParam("PrintLevel", 0);

    hephaestus::InputParameters params;
    params.SetParam("TimeStep", float(0.05));
    params.SetParam("StartTime", float(0.00));
    params.SetParam("EndTime", float(0.05));
    params.SetParam("VisualisationSteps", int(1));
    params.SetParam("UseGLVis", true);

    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("DomainProperties", domain_properties);
    params.SetParam("AuxSolvers", auxsolvers);
    params.SetParam("Postprocessors", postprocessors);
    params.SetParam("Sources", sources);
    params.SetParam("Outputs", outputs);
    params.SetParam("SolverOptions", solver_options);

    return params;
  }
};

TEST_F(TestHFormSource, CheckRun) {
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
    hephaestus::TransientFormulation *formulation =
        new hephaestus::HFormulation();
    hephaestus::TransientProblemBuilder *problem_builder =
        new hephaestus::TransientProblemBuilder(params);

    problem_builder->SetFormulation(formulation);
    problem_builder->SetMesh(pmesh);
    problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P2"));
    problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P2"));
    problem_builder->AddGridFunction(std::string("analytic_magnetic_field"),
                                     std::string("HCurl"));

    hephaestus::ProblemBuildSequencer sequencer(problem_builder);
    sequencer.ConstructEquationSystemProblem();
    std::unique_ptr<hephaestus::TransientProblem> problem =
        problem_builder->ReturnProblem();

    hephaestus::InputParameters exec_params;
    exec_params.SetParam("TimeStep", float(0.05));
    exec_params.SetParam("StartTime", float(0.00));
    exec_params.SetParam("EndTime", float(0.05));
    exec_params.SetParam("VisualisationSteps", int(1));
    exec_params.SetParam("UseGLVis", false);
    exec_params.SetParam("Problem", problem.get());
    hephaestus::TransientExecutioner *executioner =
        new hephaestus::TransientExecutioner(exec_params);

    executioner->Init();
    executioner->Execute();
    delete formulation;
  }

  hephaestus::L2ErrorVectorPostprocessor l2errpostprocessor =
      *(dynamic_cast<hephaestus::L2ErrorVectorPostprocessor *>(
          params.GetParam<hephaestus::AuxSolvers>("Postprocessors")
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
