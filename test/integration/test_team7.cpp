// Based on an H form MMS test provided by Joseph Dean

#include "auxkernels.hpp"
#include "executioner.hpp"

#include "hephaestus_transient.hpp"
#include "postprocessors.hpp"

#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestTeam7 : public testing::Test {
protected:
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
  static void source_current(const mfem::Vector &xv, double t,
                             mfem::Vector &J) {
    double x0(194e-3); // Coil centre x coordinate
    double y0(100e-3); // Coil centre y coordinate
    double a(50e-3);   // Coil thickness
    double I0(2742);   // Coil current in Ampere-turns
    double S(2.5e-3);  // Coil cross sectional area
    double freq(50.0); // Frequency in Hz

    double Jmag = (I0 / S) * sin(2 * M_PI * freq * t);

    double x = xv(0);
    double y = xv(1);
    if (abs(x - x0) < a) {
      J(0) = -(y - y0) / abs(y - y0);
    } else if (abs(y - y0) < a) {
      J(0) = 0.0;
    } else {
      J(0) = -(y - (y0 + a * ((y - y0) / abs(y - y0)))) /
             hypot(x - (x0 + a * ((x - x0) / abs(x - x0))),
                   y - (y0 + a * ((y - y0) / abs(y - y0))));
    }

    if (abs(y - y0) < a) {
      J(1) = (x - x0) / abs(x - x0);
    } else if (abs(x - x0) < a) {
      J(1) = 0.0;
    } else {
      J(1) = (x - (x0 + a * ((x - x0) / abs(x - x0)))) /
             hypot(x - (x0 + a * ((x - x0) / abs(x - x0))),
                   y - (y0 + a * ((y - y0) / abs(y - y0))));
    }

    // J(0) = 0.0;
    // // sin(M_PI * x(1)) * sin(M_PI * x(2)) * (t * 2 * M_PI * M_PI + 1);
    // J(1) = 0.0;
    J(2) = 0.0;
    J *= Jmag;
  }

  hephaestus::InputParameters test_params() {
    hephaestus::Subdomain air("air", 1);
    air.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::Subdomain plate("plate", 2);
    plate.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(3.526e7);
    hephaestus::Subdomain coil1("coil1", 3);
    coil1.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::Subdomain coil2("coil2", 4);
    coil2.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::Subdomain coil3("coil3", 5);
    coil3.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::Subdomain coil4("coil4", 6);
    coil4.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);

    hephaestus::DomainProperties domain_properties(
        std::vector<hephaestus::Subdomain>(
            {air, plate, coil1, coil2, coil3, coil4}));

    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::ConstantCoefficient(M_PI * 4.0e-7);

    hephaestus::BCMap bc_map;
    // mfem::VectorFunctionCoefficient *adotVecCoef =
    //     new mfem::VectorFunctionCoefficient(3, adot_bc);
    // bc_map["tangential_dAdt"] = new hephaestus::VectorFunctionDirichletBC(
    //     std::string("magnetic_vector_potential"), mfem::Array<int>({1, 2,
    //     3}), adotVecCoef);
    // domain_properties.vector_property_map["surface_tangential_dAdt"] =
    //     adotVecCoef;
    // domain_properties.scalar_property_map["electrical_conductivity"] =
    //     new mfem::ConstantCoefficient(1.0);

    // mfem::Array<int> ground_terminal(1);
    // ground_terminal[0] = 1;
    // mfem::FunctionCoefficient *ground_coeff =
    //     new mfem::FunctionCoefficient(potential_ground);
    // bc_map["ground_potential"] = new hephaestus::FunctionDirichletBC(
    //     std::string("electric_potential"), mfem::Array<int>({1, 2, 3}),
    //     ground_coeff);
    // domain_properties.scalar_property_map["ground_potential"] = ground_coeff;

    mfem::VectorFunctionCoefficient *JSrcCoef =
        new mfem::VectorFunctionCoefficient(3, source_current);

    mfem::Array<mfem::VectorCoefficient *> sourcecoefs(4);
    sourcecoefs[0] = JSrcCoef;
    sourcecoefs[1] = JSrcCoef;
    sourcecoefs[2] = JSrcCoef;
    sourcecoefs[3] = JSrcCoef;

    mfem::Array<int> coilsegments(4);
    coilsegments[0] = 3;
    coilsegments[1] = 4;
    coilsegments[2] = 5;
    coilsegments[3] = 6;

    mfem::PWVectorCoefficient *JSrcRestricted =
        new mfem::PWVectorCoefficient(3, coilsegments, sourcecoefs);

    domain_properties.vector_property_map["unblockedsource"] = JSrcCoef;
    domain_properties.vector_property_map["source"] = JSrcRestricted;

    mfem::VectorFunctionCoefficient *A_exact =
        new mfem::VectorFunctionCoefficient(3, A_exact_expr);
    domain_properties.vector_property_map["a_exact_coeff"] = A_exact;

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./team7_small.g")).c_str(), 1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("AVFormVisIt");
    data_collections["ParaViewDataCollection"] =
        new mfem::ParaViewDataCollection("AVFormParaView");
    hephaestus::Outputs outputs(data_collections);

    hephaestus::InputParameters hcurlvarparams;
    hcurlvarparams.SetParam("VariableName",
                            std::string("analytic_vector_potential"));
    hcurlvarparams.SetParam("FESpaceName", std::string("HCurl"));
    hcurlvarparams.SetParam("FESpaceType", std::string("Nedelec"));
    hcurlvarparams.SetParam("order", 2);
    hcurlvarparams.SetParam("components", 3);
    hephaestus::Variables variables;
    variables.AddVariable(hcurlvarparams);

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

    hephaestus::InputParameters exec_params;
    exec_params.SetParam("TimeStep", float(0.001));
    exec_params.SetParam("StartTime", float(0.00));
    exec_params.SetParam("EndTime", float(0.05));
    hephaestus::TransientExecutioner *executioner =
        new hephaestus::TransientExecutioner(exec_params);

    hephaestus::InputParameters params;
    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("Executioner", executioner);
    params.SetParam("Order", 1);
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("DomainProperties", domain_properties);
    params.SetParam("Variables", variables);
    params.SetParam("AuxKernels", auxkernels);
    params.SetParam("Postprocessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("FormulationName", std::string("AVForm"));
    std::cout << "Created params ";
    return params;
  }
};

TEST_F(TestTeam7, CheckRun) {
  hephaestus::InputParameters params(test_params());

  hephaestus::TransientExecutioner *executioner(
      params.GetParam<hephaestus::TransientExecutioner *>("Executioner"));
  std::cout << "Created exec ";
  executioner->Solve(params);

  //   hephaestus::L2ErrorVectorPostprocessor l2errpostprocessor =
  //       *(dynamic_cast<hephaestus::L2ErrorVectorPostprocessor *>(
  //           params.GetParam<hephaestus::Postprocessors>("Postprocessors")
  //               .Get("L2ErrorPostprocessor")));

  //   double r;
  //   for (std::size_t i = 1; i < l2errpostprocessor.ndofs.Size(); ++i) {
  //     r = estimate_convergence_rate(
  //         l2errpostprocessor.ndofs[i], l2errpostprocessor.ndofs[i - 1],
  //         l2errpostprocessor.l2_errs[i], l2errpostprocessor.l2_errs[i - 1],
  //         3);
  //     std::cout << r << std::endl;
  //     ASSERT_TRUE(r > params.GetParam<int>("Order"));
  //     ASSERT_TRUE(r < params.GetParam<int>("Order") + 1.0);
  //   }
}
