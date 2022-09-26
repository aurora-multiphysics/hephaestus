#include "hephaestus_transient.hpp"
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
    mfem::VectorFunctionCoefficient *hdotVecCoef =
        new mfem::VectorFunctionCoefficient(3, hdot_bc);
    bc_map["tangential_dHdt"] = new hephaestus::VectorFunctionDirichletBC(
        std::string("magnetic_field"), mfem::Array<int>({1, 2, 3}),
        hdotVecCoef);
    domain_properties.vector_property_map["surface_tangential_dHdt"] =
        hdotVecCoef;
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::ConstantCoefficient(1.0);

    mfem::Array<int> high_terminal(1);
    high_terminal[0] = 1;
    bc_map["high_potential"] = new hephaestus::FunctionDirichletBC(
        std::string("magnetic_potential"), high_terminal,
        new mfem::FunctionCoefficient(potential_high));

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 2;
    bc_map["ground_potential"] = new hephaestus::FunctionDirichletBC(
        std::string("magnetic_potential"), ground_terminal,
        new mfem::FunctionCoefficient(potential_ground));

    hephaestus::InputParameters exec_params;
    exec_params.SetParam("TimeStep", float(0.5));
    exec_params.SetParam("StartTime", float(0.0));
    exec_params.SetParam("EndTime", float(2.5));
    hephaestus::TransientExecutioner *executioner =
        new hephaestus::TransientExecutioner(exec_params);

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./cylinder-hex-q2.gen")).c_str(),
        1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("HFormVisIt");
    data_collections["ParaViewDataCollection"] =
        new mfem::ParaViewDataCollection("HFormParaView");
    hephaestus::Outputs outputs(data_collections);

    hephaestus::Variables variables;
    hephaestus::AuxKernels auxkernels;
    hephaestus::Postprocessors postprocessors;

    hephaestus::InputParameters params;
    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("Executioner", executioner);
    params.SetParam("Order", 2);
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("DomainProperties", domain_properties);
    params.SetParam("Variables", variables);
    params.SetParam("AuxKernels", auxkernels);
    params.SetParam("Postprocessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("FormulationName", std::string("HForm"));

    return params;
  }
};

TEST_F(TestHFormRod, CheckRun) {
  hephaestus::InputParameters params(test_params());
  hephaestus::TransientExecutioner *executioner(
      params.GetParam<hephaestus::TransientExecutioner *>("Executioner"));
  executioner->Init(params);
  executioner->Solve();
}
