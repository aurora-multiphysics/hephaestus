#include "hephaestus_transient.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class CopperConductivityCoefficient : public hephaestus::CoupledCoefficient {
public:
  CopperConductivityCoefficient(const hephaestus::InputParameters &params)
      : hephaestus::CoupledCoefficient(params){};
  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip) {
    return 2.0 * M_PI * 10 + 0.001 * (gf->GetValue(T, ip));
  }
};

class TestEBFormCoupled : public testing::Test {
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

    // materialCopper instances and get property coefs? init can be for all...
    // CoupledCoefficients must also be added to Postprocessors
    hephaestus::InputParameters copper_conductivity_params;
    copper_conductivity_params.SetParam("CoupledVariableName",
                                        std::string("electric_potential"));

    hephaestus::CoupledCoefficient *wireConductivity =
        new CopperConductivityCoefficient(copper_conductivity_params);

    hephaestus::Postprocessors postprocessors;
    postprocessors.Register("CoupledCoefficient", wireConductivity, false);

    hephaestus::Subdomain wire("wire", 1);
    wire.property_map["electrical_conductivity"] = wireConductivity;

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
    bc_map["tangential_dEdt"] = new hephaestus::VectorFunctionDirichletBC(
        std::string("electric_field"), mfem::Array<int>({1, 2, 3}),
        edotVecCoef);
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::ConstantCoefficient(1.0);
    domain_properties.vector_property_map["surface_tangential_dEdt"] =
        edotVecCoef;

    mfem::Array<int> high_terminal(1);
    high_terminal[0] = 1;
    mfem::FunctionCoefficient *potential_src =
        new mfem::FunctionCoefficient(potential_high);
    bc_map["high_potential"] = new hephaestus::FunctionDirichletBC(
        std::string("electric_potential"), high_terminal, potential_src);
    domain_properties.scalar_property_map["source_potential"] = potential_src;

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 2;
    bc_map["ground_potential"] = new hephaestus::FunctionDirichletBC(
        std::string("electric_potential"), ground_terminal,
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
        new mfem::VisItDataCollection("EBFormVisIt");
    data_collections["ParaViewDataCollection"] =
        new mfem::ParaViewDataCollection("EBFormParaView");
    hephaestus::Outputs outputs(data_collections);

    hephaestus::Variables variables;
    hephaestus::AuxKernels auxkernels;

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
    params.SetParam("FormulationName", std::string("EBForm"));

    return params;
  }
};

TEST_F(TestEBFormCoupled, CheckRun) {
  hephaestus::InputParameters params(test_params());
  hephaestus::TransientExecutioner *executioner(
      params.GetParam<hephaestus::TransientExecutioner *>("Executioner"));
  executioner->Init(params);
  executioner->Solve();
}
