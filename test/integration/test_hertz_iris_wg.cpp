#include "hephaestus_hertz.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestHertzIrisWaveguide : public testing::Test {
protected:
  static void e_bc_r(const mfem::Vector &x, mfem::Vector &E) {
    E.SetSize(3);
    E = 0.0;
  }

  static void e_bc_i(const mfem::Vector &x, mfem::Vector &E) {
    E.SetSize(3);
    E = 0.0;
  }

  hephaestus::Inputs hertz_example_inputs() {
    hephaestus::BCMap bc_map;

    hephaestus::VectorFunctionDirichletBC e_bc(
        std::string("boundary_1"), mfem::Array<int>({1, 2}),
        new mfem::VectorFunctionCoefficient(3, e_bc_r),
        new mfem::VectorFunctionCoefficient(3, e_bc_i));

    bc_map["tangential_E"] = new hephaestus::VectorFunctionDirichletBC(e_bc);

    hephaestus::Subdomain air("air", 1);

    air.property_map["real_electrical_conductivity"] =
        new mfem::ConstantCoefficient(0.0);
    air.property_map["imag_electrical_conductivity"] =
        new mfem::ConstantCoefficient(0.0);
    air.property_map["real_rel_permittivity"] =
        new mfem::ConstantCoefficient(1.0);
    air.property_map["imag_rel_permittivity"] =
        new mfem::ConstantCoefficient(0.0);
    air.property_map["real_rel_permeability"] =
        new mfem::ConstantCoefficient(1.0);
    air.property_map["imag_rel_permeability"] =
        new mfem::ConstantCoefficient(0.0);

    hephaestus::DomainProperties material_map(
        std::vector<hephaestus::Subdomain>({air}));

    hephaestus::Executioner executioner(std::string("transient"), 0.5, 5.0);
    mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./irises.g")).c_str(),
                    1, 1);


    std::map<std::string, mfem::DataCollection*> data_collections;
    data_collections["VisItDataCollection"] = new mfem::VisItDataCollection("Hertz-AMR-Parallel");
    hephaestus::Outputs outputs(data_collections);
    hephaestus::Inputs inputs(mesh, std::string("Hertz"), 2, bc_map,
                              material_map, executioner, outputs);
    return inputs;
  }
};

TEST_F(TestHertzIrisWaveguide, CheckRun) {
  hephaestus::Inputs inputs(hertz_example_inputs());
  std::vector<char *> argv;
  hertz_solve(0, argv.data(), inputs);
}
