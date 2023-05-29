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
    wire.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::Subdomain air("air", 2);
    air.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);

    hephaestus::DomainProperties domain_properties(
        std::vector<hephaestus::Subdomain>({wire, air}));

    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::FunctionCoefficient(mu_expr);

    hephaestus::BCMap bc_map;
    mfem::VectorFunctionCoefficient *adotVecCoef =
        new mfem::VectorFunctionCoefficient(3, adot_bc);
    bc_map.Register("tangential_dAdt",
                    new hephaestus::VectorFunctionDirichletBC(
                        std::string("dmagnetic_vector_potential_dt"),
                        mfem::Array<int>({1, 2, 3}), adotVecCoef),
                    true);
    domain_properties.vector_property_map["surface_tangential_dAdt"] =
        adotVecCoef;
    domain_properties.scalar_property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 1;
    mfem::FunctionCoefficient *ground_coeff =
        new mfem::FunctionCoefficient(potential_ground);
    bc_map.Register("ground_potential",
                    new hephaestus::FunctionDirichletBC(
                        std::string("electric_potential"),
                        mfem::Array<int>({1, 2, 3}), ground_coeff),
                    true);
    domain_properties.scalar_property_map["ground_potential"] = ground_coeff;

    mfem::VectorFunctionCoefficient *A_exact =
        new mfem::VectorFunctionCoefficient(3, A_exact_expr);
    domain_properties.vector_property_map["a_exact_coeff"] = A_exact;

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./beam-tet.mesh")).c_str(), 1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("AVFormVisIt");
    hephaestus::Outputs outputs(data_collections);

    hephaestus::InputParameters hcurlfespaceparams;
    hcurlfespaceparams.SetParam("FESpaceName", std::string("HCurl"));
    hcurlfespaceparams.SetParam("FESpaceType", std::string("ND"));
    hcurlfespaceparams.SetParam("order", var_order);
    hcurlfespaceparams.SetParam("components", 3);
    hephaestus::FESpaces fespaces;
    fespaces.StoreInput(hcurlfespaceparams);

    hephaestus::InputParameters analyicaparams;
    analyicaparams.SetParam("VariableName",
                            std::string("analytic_vector_potential"));
    analyicaparams.SetParam("FESpaceName", std::string("HCurl"));
    hephaestus::GridFunctions gridfunctions;
    gridfunctions.StoreInput(analyicaparams);

    hephaestus::InputParameters l2errpostprocparams;
    l2errpostprocparams.SetParam("VariableName",
                                 std::string("magnetic_vector_potential"));
    l2errpostprocparams.SetParam("VectorCoefficientName",
                                 std::string("a_exact_coeff"));
    hephaestus::Postprocessors postprocessors;
    postprocessors.Register(
        "L2ErrorPostprocessor",
        new hephaestus::L2ErrorVectorPostprocessor(l2errpostprocparams), true);

    hephaestus::InputParameters vectorcoeffauxparams;
    vectorcoeffauxparams.SetParam("VariableName",
                                  std::string("analytic_vector_potential"));
    vectorcoeffauxparams.SetParam("VectorCoefficientName",
                                  std::string("a_exact_coeff"));

    hephaestus::AuxKernels auxkernels;
    auxkernels.Register(
        "VectorCoefficientAuxKernel",
        new hephaestus::VectorCoefficientAuxKernel(vectorcoeffauxparams), true);

    hephaestus::Sources sources;
    mfem::VectorFunctionCoefficient *JSrcCoef =
        new mfem::VectorFunctionCoefficient(3, source_field);
    domain_properties.vector_property_map["source"] = JSrcCoef;
    hephaestus::InputParameters div_free_source_params;
    div_free_source_params.SetParam("SourceName", std::string("source"));
    div_free_source_params.SetParam("HCurlFESpaceName",
                                    std::string("_HCurlFESpace"));
    div_free_source_params.SetParam("H1FESpaceName", std::string("_H1FESpace"));
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
    params.SetParam("DomainProperties", domain_properties);
    params.SetParam("FESpaces", fespaces);
    params.SetParam("GridFunctions", gridfunctions);
    params.SetParam("AuxKernels", auxkernels);
    params.SetParam("Postprocessors", postprocessors);
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

    mfem::ParMesh pmesh(unrefined_pmesh);

    for (int l = 0; l < par_ref_levels; l++) {
      pmesh.UniformRefinement();
    }
    params.SetParam("Mesh", pmesh);
    hephaestus::TransientFormulation *formulation =
        new hephaestus::AVFormulation();
    params.SetParam("Formulation", formulation);
    hephaestus::TransientProblemBuilder *problem_builder =
        new hephaestus::TransientProblemBuilder(params);
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
  }

  hephaestus::L2ErrorVectorPostprocessor l2errpostprocessor =
      *(dynamic_cast<hephaestus::L2ErrorVectorPostprocessor *>(
          params.GetParam<hephaestus::Postprocessors>("Postprocessors")
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
