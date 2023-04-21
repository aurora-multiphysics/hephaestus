#include "hephaestus_joule.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestJouleRod : public testing::Test {
protected:
  static double potential(const mfem::Vector &x, double t) {
    double wj_(2.0 * M_PI / 60.0);
    // the value
    double T;
    if (x[2] < 0.0) {
      T = 1.0;
    } else {
      T = -1.0;
    }

    return T * cos(wj_ * t);
  }

  hephaestus::Inputs joule_rod_inputs() {
    hephaestus::BCMap bc_map;
    bc_map.Register("tangential_dEdt",
                    new hephaestus::BoundaryCondition(
                        std::string("boundary_1"), mfem::Array<int>({1, 2, 3})),
                    true);

    bc_map.Register("thermal_flux",
                    new hephaestus::BoundaryCondition(std::string("boundary_2"),
                                                      mfem::Array<int>({1, 2})),
                    true);

    bc_map.Register("electric_potential",
                    new hephaestus::FunctionDirichletBC(
                        std::string("boundary_3"), mfem::Array<int>({1, 2}),
                        new mfem::FunctionCoefficient(potential)),
                    true);

    double sigma = 2.0 * M_PI * 10;
    double Tcapacity = 1.0;
    double Tconductivity = 0.01;

    double sigmaAir;
    double TcondAir;
    double TcapAir;

    sigmaAir = 1.0e-6 * sigma;
    TcondAir = 1.0e6 * Tconductivity;
    TcapAir = 1.0 * Tcapacity;

    hephaestus::Subdomain wire("wire", 1);
    wire.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(sigma);
    wire.property_map["heat_capacity"] =
        new mfem::ConstantCoefficient(Tcapacity);
    wire.property_map["inverse_heat_capacity"] =
        new mfem::ConstantCoefficient(1.0 / Tcapacity);
    wire.property_map["inverse_thermal_conductivity"] =
        new mfem::ConstantCoefficient(1.0 / Tconductivity);

    hephaestus::Subdomain air("air", 2);
    air.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(sigmaAir);
    air.property_map["heat_capacity"] = new mfem::ConstantCoefficient(TcapAir);
    air.property_map["inverse_heat_capacity"] =
        new mfem::ConstantCoefficient(1.0 / TcapAir);
    air.property_map["inverse_thermal_conductivity"] =
        new mfem::ConstantCoefficient(1.0 / TcondAir);

    hephaestus::DomainProperties material_map(
        std::vector<hephaestus::Subdomain>({wire, air}));

    hephaestus::Executioner executioner(std::string("transient"), 0.5, 0.0,
                                        2.5);
    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./cylinder-hex-q2.gen")).c_str(),
        1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("JouleRodVisIt");
    data_collections["ParaViewDataCollection"] =
        new mfem::ParaViewDataCollection("JouleRodParaView");
    hephaestus::Outputs outputs(data_collections);
    hephaestus::Inputs inputs(mesh, std::string("Joule"), 2, bc_map,
                              material_map, executioner, outputs);
    return inputs;
  }
};

TEST_F(TestJouleRod, CheckRun) {
  hephaestus::Inputs inputs(joule_rod_inputs());
  std::vector<char *> argv;
  joule_solve(0, argv.data(), inputs);
}
