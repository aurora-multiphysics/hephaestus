#include "hephaestus.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestEBFormRod : public testing::Test {
protected:
  static double potential_high(const mfem::Vector &x, double t) {
    double wj_(2.0 * M_PI / 60.0);
    return 2 * cos(wj_ * t);
  }
  static double potential_ground(const mfem::Vector &x, double t) {
    return 0.0;
  }
  static void edot_bc(const mfem::Vector &x, mfem::Vector &E) { E = 0.0; }
  static void j_src(const mfem::Vector &x, double t, mfem::Vector &j) {
    double wj_(2.0 * M_PI / 60.0);
    j[0] = 0.0;
    j[1] = 0.0;
    j[2] = 2 * sin(wj_ * t);
  }

  hephaestus::InputParameters test_params() {
    double sigma = 2.0 * M_PI * 10;

    double sigmaAir;

    sigmaAir = 1.0e-6 * sigma;

    hephaestus::Subdomain wire("wire", 1);
    wire.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(sigma);

    hephaestus::Subdomain air("air", 2);
    air.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(sigmaAir);

    hephaestus::DomainProperties domain_properties(
        std::vector<hephaestus::Subdomain>({wire, air}));

    domain_properties.scalar_property_map["electrical_conductivity"] =
        new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
            std::string("electrical_conductivity")));

    hephaestus::BCMap bc_map;
    mfem::VectorFunctionCoefficient *edotVecCoef =
        new mfem::VectorFunctionCoefficient(3, edot_bc);
    bc_map.Register("tangential_dEdt",
                    new hephaestus::VectorFunctionDirichletBC(
                        std::string("electric_field"),
                        mfem::Array<int>({1, 2, 3}), edotVecCoef),
                    true);
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::ConstantCoefficient(1.0);
    domain_properties.vector_property_map["surface_tangential_dEdt"] =
        edotVecCoef;

    mfem::Array<int> high_terminal(1);
    high_terminal[0] = 1;
    mfem::FunctionCoefficient *potential_src =
        new mfem::FunctionCoefficient(potential_high);
    bc_map.Register(
        "high_potential",
        new hephaestus::FunctionDirichletBC(std::string("electric_potential"),
                                            high_terminal, potential_src),
        true);
    domain_properties.scalar_property_map["source_potential"] = potential_src;

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 2;
    bc_map.Register("ground_potential",
                    new hephaestus::FunctionDirichletBC(
                        std::string("electric_potential"), ground_terminal,
                        new mfem::FunctionCoefficient(potential_ground)),
                    true);

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./cylinder-hex-q2.gen")).c_str(),
        1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("EBFormVisIt");
    data_collections["ParaViewDataCollection"] =
        new mfem::ParaViewDataCollection("EBFormParaView");
    hephaestus::Outputs outputs(data_collections);

    hephaestus::InputParameters hcurlfespaceparams;
    hcurlfespaceparams.SetParam("FESpaceName", std::string("HCurl"));
    hcurlfespaceparams.SetParam("FESpaceType", std::string("ND"));
    hcurlfespaceparams.SetParam("order", 2);
    hcurlfespaceparams.SetParam("components", 3);
    hephaestus::InputParameters h1fespaceparams;
    h1fespaceparams.SetParam("FESpaceName", std::string("H1"));
    h1fespaceparams.SetParam("FESpaceType", std::string("H1"));
    h1fespaceparams.SetParam("order", 2);
    h1fespaceparams.SetParam("components", 3);
    hephaestus::FESpaces fespaces;
    fespaces.StoreInput(hcurlfespaceparams);
    fespaces.StoreInput(h1fespaceparams);

    hephaestus::GridFunctions gridfunctions;
    hephaestus::AuxKernels auxkernels;
    hephaestus::Postprocessors postprocessors;
    hephaestus::Sources sources;
    hephaestus::InputParameters scalar_potential_source_params;
    scalar_potential_source_params.SetParam("SourceName",
                                            std::string("source"));
    scalar_potential_source_params.SetParam("PotentialName",
                                            std::string("electric_potential"));
    scalar_potential_source_params.SetParam("HCurlFESpaceName",
                                            std::string("HCurl"));
    scalar_potential_source_params.SetParam("H1FESpaceName", std::string("H1"));
    scalar_potential_source_params.SetParam(
        "ConductivityCoefName", std::string("electrical_conductivity"));
    hephaestus::InputParameters current_solver_options;
    current_solver_options.SetParam("Tolerance", float(1.0e-9));
    current_solver_options.SetParam("MaxIter", (unsigned int)1000);
    current_solver_options.SetParam("PrintLevel", -1);
    scalar_potential_source_params.SetParam("SolverOptions",
                                            current_solver_options);
    sources.Register(
        "source",
        new hephaestus::ScalarPotentialSource(scalar_potential_source_params),
        true);

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-9));
    solver_options.SetParam("MaxIter", (unsigned int)1000);
    solver_options.SetParam("PrintLevel", 0);

    hephaestus::TransientFormulation *formulation =
        new hephaestus::EBDualFormulation();

    hephaestus::InputParameters params;
    params.SetParam("TimeStep", float(0.5));
    params.SetParam("StartTime", float(0.0));
    params.SetParam("EndTime", float(2.5));
    params.SetParam("VisualisationSteps", int(1));
    params.SetParam("UseGLVis", true);

    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("DomainProperties", domain_properties);
    params.SetParam("FESpaces", fespaces);
    params.SetParam("GridFunctions", gridfunctions);
    params.SetParam("AuxKernels", auxkernels);
    params.SetParam("Postprocessors", postprocessors);
    params.SetParam("Sources", sources);
    params.SetParam("Outputs", outputs);
    params.SetParam("Formulation", formulation);
    params.SetParam("SolverOptions", solver_options);

    return params;
  }
};

TEST_F(TestEBFormRod, CheckRun) {
  hephaestus::InputParameters params(test_params());
  hephaestus::TransientExecutioner *executioner =
      new hephaestus::TransientExecutioner(params);
  executioner->Init();
  executioner->Solve();
}
