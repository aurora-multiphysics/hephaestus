// Based on an H form MMS test provided by Joseph Dean

#include "auxkernels.hpp"
#include "executioner.hpp"

#include "hephaestus_transient.hpp"
#include "postprocessors.hpp"

#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestCurrentRing : public testing::Test {
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

  static void source_current(const mfem::Vector &x, double t, mfem::Vector &j) {

    j.SetSize(x.Size());
    j = 0.0;

    mfem::Vector a(x.Size()); // Normalized Axis vector
    a(0) = 0.5;
    a(1) = 0.5;
    a(2) = 0.45;

    mfem::Vector xu(x.Size()); // x vector relative to the axis end-point
    xu(0) = 0.5;
    xu(1) = 0.5;
    xu(2) = 0.55;

    xu *= -1.0;
    a += xu;
    xu += x;

    mfem::Vector ju(x.Size()); // Unit vector in direction of current

    double h = a.Norml2();

    if (h == 0.0) {
      return;
    }

    double ra = 0.2;
    double rb = 0.3;
    if (ra > rb) {
      double rc = ra;
      ra = rb;
      rb = rc;
    }
    double xa = xu * a;

    if (h > 0.0) {
      xu.Add(-xa / (h * h), a);
    }

    double xp = xu.Norml2();

    if (xa >= 0.0 && xa <= h * h && xp >= ra && xp <= rb) {
      ju(0) = a(1) * xu(2) - a(2) * xu(1);
      ju(1) = a(2) * xu(0) - a(0) * xu(2);
      ju(2) = a(0) * xu(1) - a(1) * xu(0);
      ju /= h;

      j.Add(1.0 / (h * (rb - ra)), ju);
    }
  }

  hephaestus::InputParameters test_params() {
    hephaestus::Subdomain coil("coil", 1);
    coil.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::Subdomain ring("ring", 3);
    ring.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(2.273e7);
    hephaestus::Subdomain air("air", 2);
    air.property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);

    hephaestus::DomainProperties domain_properties(
        std::vector<hephaestus::Subdomain>({coil, air, ring}));

    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::ConstantCoefficient(M_PI * 4.0e-7);
    domain_properties.scalar_property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(1.0);
    hephaestus::BCMap bc_map;

    mfem::VectorFunctionCoefficient *JSrcCoef =
        new mfem::VectorFunctionCoefficient(3, source_current);

    mfem::Array<int> source_coil(1);
    source_coil[0] = 1;
    mfem::Array<mfem::VectorCoefficient *> sources(1);
    sources[0] = JSrcCoef;

    mfem::PWVectorCoefficient *JRestricted =
        new mfem::PWVectorCoefficient(3, source_coil, sources);

    // mfem::VectorRestrictedCoefficient *JRestricted =
    //     new mfem::VectorRestrictedCoefficient(*JSrcCoef, source_coil);
    domain_properties.vector_property_map["source"] = JRestricted;

    mfem::VectorFunctionCoefficient *A_exact =
        new mfem::VectorFunctionCoefficient(3, A_exact_expr);

    domain_properties.vector_property_map["a_exact_coeff"] = A_exact;

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./solenoid.g")).c_str(), 1, 1);
    mesh.UniformRefinement();
    mesh.SetCurvature(2);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("AVFormRingVisIt");
    data_collections["ParaViewDataCollection"] =
        new mfem::ParaViewDataCollection("AVFormRingParaView");

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
    exec_params.SetParam("TimeStep", float(0.1));
    exec_params.SetParam("StartTime", float(0.00));
    exec_params.SetParam("EndTime", float(1.0));
    hephaestus::TransientExecutioner *executioner =
        new hephaestus::TransientExecutioner(exec_params);

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
    params.SetParam("FormulationName", std::string("AVForm"));

    return params;
  }
};

TEST_F(TestCurrentRing, CheckRun) {
  hephaestus::InputParameters params(test_params());

  hephaestus::TransientExecutioner *executioner(
      params.GetParam<hephaestus::TransientExecutioner *>("Executioner"));
  executioner->Solve(params);
}
