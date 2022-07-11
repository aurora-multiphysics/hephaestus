#include "hephaestus_transient.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestAVForm : public testing::Test {
protected:
  static double potential_high(const mfem::Vector &x, double t) {
    double wj_(2.0 * M_PI / 60.0);
    return 2 * cos(wj_ * t);
  }
  static double potential_ground(const mfem::Vector &x, double t) {
    return 0.0;
  }
  static void adot_bc(const mfem::Vector &x, mfem::Vector &E) { E = 0.0; }
  static void j_src(const mfem::Vector &x, double t, mfem::Vector &j) {
    double wj_(2.0 * M_PI / 60.0);
    j[0] = 0.0;
    j[1] = 0.0;
    j[2] = 2 * sin(wj_ * t);
  }

  hephaestus::Inputs avform_rod_inputs() {
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
    mfem::VectorFunctionCoefficient *adotVecCoef =
        new mfem::VectorFunctionCoefficient(3, adot_bc);
    bc_map["tangential_dAdt"] = new hephaestus::VectorFunctionDirichletBC(
        std::string("magnetic_vector_potential"), mfem::Array<int>({1, 2, 3}),
        adotVecCoef);
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::ConstantCoefficient(1.0);
    domain_properties.vector_property_map["surface_tangential_dAdt"] =
        adotVecCoef;

    mfem::Array<int> high_terminal(1);
    high_terminal[0] = 1;
    // mfem::VectorFunctionCoefficient *jVecCoef =
    //     new mfem::VectorFunctionCoefficient(3, j_src);
    // bc_map["current_inlet"] = new hephaestus::IntegratedBC(
    //     std::string("electric_potential"), high_terminal,
    //     new mfem::BoundaryNormalLFIntegrator(*jVecCoef));
    // domain_properties.vector_property_map["surface_normal_current_density"] =
    //     jVecCoef;
    bc_map["high_potential"] = new hephaestus::FunctionDirichletBC(
        std::string("electric_potential"), high_terminal,
        new mfem::FunctionCoefficient(potential_high));

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 2;
    bc_map["ground_potential"] = new hephaestus::FunctionDirichletBC(
        std::string("electric_potential"), ground_terminal,
        new mfem::FunctionCoefficient(potential_ground));

    hephaestus::Executioner executioner(std::string("transient"), 0.5, 0.0,
                                        2.5);
    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./cylinder-hex-q2.gen")).c_str(),
        1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("AVFormVisIt");
    data_collections["ParaViewDataCollection"] =
        new mfem::ParaViewDataCollection("AVFormParaview");
    hephaestus::Outputs outputs(data_collections);
    hephaestus::Inputs inputs(mesh, std::string("AVForm"), 2, bc_map,
                              domain_properties, executioner, outputs);
    return inputs;
  }
};

TEST_F(TestAVForm, CheckRun) {
  hephaestus::Inputs inputs(avform_rod_inputs());
  std::vector<char *> argv;
  transient_solve(0, argv.data(), inputs);
}
