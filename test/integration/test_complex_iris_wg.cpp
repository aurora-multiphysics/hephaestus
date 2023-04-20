#include "../common/pfem_extras.hpp"
#include "factory.hpp"
#include "inputs.hpp"
#include "steady_executioner.hpp"
#include <gtest/gtest.h>

extern const char *DATA_DIR;

class TestComplexIrisWaveguide : public testing::Test {
protected:
  static void zero_func(const mfem::Vector &x, mfem::Vector &v) {
    v.SetSize(3);
    v = 0.0;
  }

  static void e_bc_r(const mfem::Vector &x, mfem::Vector &E) {
    E.SetSize(3);
    E = 0.0;
  }

  static void e_bc_i(const mfem::Vector &x, mfem::Vector &E) {
    E.SetSize(3);
    E = 0.0;
  }

  // inline static const double epsilon0_ = 8.8541878176e-12;
  inline static const double epsilon0_ = 8.8541878176e-12;
  // Permeability of Free Space (units H/m)
  inline static const double mu0_ = 4.0e-7 * M_PI;
  inline static const double freq_ = 9.3e9; // 10/2pi
  double port_length_vector[3] = {0.0, 22.86e-3, 0.0};
  double port_width_vector[3] = {0.0, 0.0, 10.16e-3};

  mfem::Array<int> dirichlet_attr;

  static void RWTE10(const mfem::Vector &x,
                     std::vector<std::complex<double>> &E) {
    std::complex<double> zi = std::complex<double>(0., 1.);
    double port_length_vector[3] = {0.0, 22.86e-3, 0.0};
    double port_width_vector[3] = {0.0, 0.0, 10.16e-3};
    double omega_ = 2 * M_PI * freq_;
    mfem::Vector a1Vec(port_length_vector, 3);
    mfem::Vector a2Vec(port_width_vector, 3);
    mfem::Vector a3Vec;
    mfem::Vector a2xa3;
    mfem::Vector a3xa1;

    hephaestus::cross_product(a1Vec, a2Vec, a3Vec);
    hephaestus::cross_product(a2Vec, a3Vec, a2xa3);
    hephaestus::cross_product(a3Vec, a1Vec, a3xa1);

    mfem::Vector k_a = a2xa3;
    mfem::Vector k_c = a3Vec;
    double V = InnerProduct(a1Vec, a2xa3);

    double kc = M_PI / a1Vec.Norml2();
    double k0 = omega_ * sqrt(epsilon0_ * mu0_);
    std::complex<double> k_ = std::complex<double>(0., sqrt(k0 * k0 - kc * kc));
    k_a *= M_PI / V;
    k_c *= k_.imag() / a3Vec.Norml2();

    mfem::Vector E_hat;
    hephaestus::cross_product(k_c, k_a, E_hat);
    E_hat *= 1.0 / E_hat.Norml2();

    double E0(sqrt(2 * omega_ * mu0_ /
                   (a1Vec.Norml2() * a2Vec.Norml2() * k_.imag())));
    std::complex<double> E_mag =
        E0 * sin(InnerProduct(k_a, x)) * exp(-zi * InnerProduct(k_c, x));

    E[0] = E_mag * E_hat(1);
    E[1] = E_mag * E_hat(2);
    E[2] = E_mag * E_hat(0);
  }

  static void RWTE10_real(const mfem::Vector &x, mfem::Vector &v) {
    std::vector<std::complex<double>> Eval(x.Size());

    std::complex<double> zi = std::complex<double>(0., 1.);
    double port_length_vector[3] = {0.0, 22.86e-3, 0.0};
    double port_width_vector[3] = {0.0, 0.0, 10.16e-3};
    double omega_ = 2 * M_PI * freq_;
    mfem::Vector a1Vec(port_length_vector, 3);
    mfem::Vector a2Vec(port_width_vector, 3);
    mfem::Vector a3Vec;
    mfem::Vector a2xa3;
    mfem::Vector a3xa1;

    hephaestus::cross_product(a1Vec, a2Vec, a3Vec);
    hephaestus::cross_product(a2Vec, a3Vec, a2xa3);
    hephaestus::cross_product(a3Vec, a1Vec, a3xa1);

    mfem::Vector k_a = a2xa3;
    mfem::Vector k_c = a3Vec;
    double V = InnerProduct(a1Vec, a2xa3);

    double kc = M_PI / a1Vec.Norml2();
    double k0 = omega_ * sqrt(epsilon0_ * mu0_);
    std::complex<double> k_ = std::complex<double>(0., sqrt(k0 * k0 - kc * kc));

    RWTE10(x, Eval);
    for (int i = 0; i < x.Size(); ++i) {
      v(i) = -2 * k_.imag() * Eval[i].imag() / mu0_;
    }
  }
  static void RWTE10_imag(const mfem::Vector &x, mfem::Vector &v) {
    std::vector<std::complex<double>> Eval(x.Size());

    std::complex<double> zi = std::complex<double>(0., 1.);
    double port_length_vector[3] = {0.0, 22.86e-3, 0.0};
    double port_width_vector[3] = {0.0, 0.0, 10.16e-3};
    double omega_ = 2 * M_PI * freq_;
    mfem::Vector a1Vec(port_length_vector, 3);
    mfem::Vector a2Vec(port_width_vector, 3);
    mfem::Vector a3Vec;
    mfem::Vector a2xa3;
    mfem::Vector a3xa1;

    hephaestus::cross_product(a1Vec, a2Vec, a3Vec);
    hephaestus::cross_product(a2Vec, a3Vec, a2xa3);
    hephaestus::cross_product(a3Vec, a1Vec, a3xa1);

    mfem::Vector k_a = a2xa3;
    mfem::Vector k_c = a3Vec;
    double V = InnerProduct(a1Vec, a2xa3);

    double kc = M_PI / a1Vec.Norml2();
    double k0 = omega_ * sqrt(epsilon0_ * mu0_);
    std::complex<double> k_ = std::complex<double>(0., sqrt(k0 * k0 - kc * kc));

    RWTE10(x, Eval);
    for (int i = 0; i < x.Size(); ++i) {
      v(i) = 2 * k_.imag() * Eval[i].real() / mu0_;
    }
  }

  hephaestus::InputParameters test_params() {
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

    hephaestus::DomainProperties domain_properties(
        std::vector<hephaestus::Subdomain>({air}));

    domain_properties.scalar_property_map["frequency"] =
        new mfem::ConstantCoefficient(9.3e9);
    domain_properties.vector_property_map["UReal"] =
        new mfem::VectorFunctionCoefficient(3, RWTE10_real);
    domain_properties.vector_property_map["UImag"] =
        new mfem::VectorFunctionCoefficient(3, RWTE10_imag);
    domain_properties.vector_property_map["Zero"] =
        new mfem::VectorFunctionCoefficient(3, zero_func);

    hephaestus::BCMap bc_map;
    // mfem::Array<int> dirichlet_attr(1);
    dirichlet_attr.Append(1);
    bc_map["tangential_E"] = new hephaestus::VectorFunctionDirichletBC(
        std::string("electric_field"), dirichlet_attr,
        new mfem::VectorFunctionCoefficient(3, e_bc_r),
        new mfem::VectorFunctionCoefficient(3, e_bc_i));

    mfem::Vector a1Vec(port_length_vector, 3);
    double kc = M_PI / a1Vec.Norml2();
    double omega_ = 2.0 * M_PI * freq_;
    double k0 = omega_ * sqrt(epsilon0_ * mu0_);
    std::complex<double> k_ = std::complex<double>(0., sqrt(k0 * k0 - kc * kc));
    // Robin coefficient, but already handled here.
    // domain_properties.scalar_property_map["etaInv"] =
    //     new mfem::ConstantCoefficient(-k_.imag() / (mu0_ * omega_));
    domain_properties.scalar_property_map["abc"] =
        new mfem::ConstantCoefficient(k_.imag() / mu0_);
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::ConstantCoefficient(mu0_);
    domain_properties.scalar_property_map["dielectric_permittivity"] =
        new mfem::ConstantCoefficient(epsilon0_);
    domain_properties.scalar_property_map["electrical_conductivity"] =
        new mfem::ConstantCoefficient(0.0);

    mfem::Array<int> wgi_in_attr(1);
    wgi_in_attr[0] = 2;
    bc_map["WaveguidePortIn"] = new hephaestus::RobinBC(
        std::string("electric_field"), wgi_in_attr, NULL,
        new mfem::VectorFEBoundaryTangentLFIntegrator(
            *(domain_properties.vector_property_map["UReal"])),
        new mfem::VectorFEMassIntegrator(
            *domain_properties.scalar_property_map["abc"]),
        new mfem::VectorFEBoundaryTangentLFIntegrator(
            *(domain_properties.vector_property_map["UImag"])));

    mfem::Array<int> wgi_out_attr(1);
    wgi_out_attr[0] = 3;
    bc_map["WaveguidePortOut"] = new hephaestus::RobinBC(
        std::string("electric_field"), wgi_out_attr, NULL,
        new mfem::VectorFEBoundaryTangentLFIntegrator(
            *(domain_properties.vector_property_map["Zero"])),
        new mfem::VectorFEMassIntegrator(
            *domain_properties.scalar_property_map["abc"]),
        new mfem::VectorFEBoundaryTangentLFIntegrator(
            *(domain_properties.vector_property_map["Zero"])));

    mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./irises.g")).c_str(),
                    1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("Hertz-AMR-Parallel-VisIt");
    hephaestus::Outputs outputs(data_collections);

    // hephaestus::InputParameters hcurlfespaceparams;
    // hcurlfespaceparams.SetParam("FESpaceName", std::string("HCurl"));
    // hcurlfespaceparams.SetParam("FESpaceType", std::string("ND"));
    // hcurlfespaceparams.SetParam("order", 2);
    // hcurlfespaceparams.SetParam("components", 3);
    // hephaestus::InputParameters h1fespaceparams;
    // h1fespaceparams.SetParam("FESpaceName", std::string("H1"));
    // h1fespaceparams.SetParam("FESpaceType", std::string("H1"));
    // h1fespaceparams.SetParam("order", 2);
    // h1fespaceparams.SetParam("components", 3);
    hephaestus::FESpaces fespaces;
    // fespaces.StoreInput(hcurlfespaceparams);
    // fespaces.StoreInput(h1fespaceparams);

    // hephaestus::InputParameters analyicaparams;
    // analyicaparams.SetParam("VariableName",
    //                         std::string("analytic_vector_potential"));
    // analyicaparams.SetParam("FESpaceName", std::string("HCurl"));
    hephaestus::GridFunctions gridfunctions;
    // gridfunctions.StoreInput(analyicaparams);

    // hephaestus::InputParameters l2errpostprocparams;
    // l2errpostprocparams.SetParam("VariableName",
    //                              std::string("magnetic_vector_potential"));
    // l2errpostprocparams.SetParam("VectorCoefficientName",
    //                              std::string("a_exact_coeff"));
    hephaestus::Postprocessors postprocessors;
    // postprocessors.Register(
    //     "L2ErrorPostprocessor",
    //     new hephaestus::L2ErrorVectorPostprocessor(l2errpostprocparams),
    //     true);

    // hephaestus::InputParameters vectorcoeffauxparams;
    // vectorcoeffauxparams.SetParam("VariableName",
    //                               std::string("analytic_vector_potential"));
    // vectorcoeffauxparams.SetParam("VectorCoefficientName",
    //                               std::string("a_exact_coeff"));

    hephaestus::AuxKernels auxkernels;
    // auxkernels.Register(
    //     "VectorCoefficientAuxKernel",
    //     new hephaestus::VectorCoefficientAuxKernel(vectorcoeffauxparams),
    //     true);

    hephaestus::Sources sources;
    // mfem::VectorFunctionCoefficient *JSrcCoef =
    //     new mfem::VectorFunctionCoefficient(3, source_field);
    // domain_properties.vector_property_map["source"] = JSrcCoef;
    // hephaestus::InputParameters div_free_source_params;
    // div_free_source_params.SetParam("SourceName", std::string("source"));
    // div_free_source_params.SetParam("HCurlFESpaceName",
    //                                 std::string("_HCurlFESpace"));
    // div_free_source_params.SetParam("H1FESpaceName", std::string("H1"));
    // hephaestus::InputParameters current_solver_options;
    // current_solver_options.SetParam("Tolerance", float(1.0e-12));
    // current_solver_options.SetParam("MaxIter", (unsigned int)200);
    // current_solver_options.SetParam("PrintLevel", 0);
    // div_free_source_params.SetParam("SolverOptions", current_solver_options);
    // sources.Register(
    //     "source",
    //     new hephaestus::DivFreeVolumetricSource(div_free_source_params),
    //     true);

    hephaestus::HertzFormulation *formulation =
        new hephaestus::HertzFormulation();

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-16));
    solver_options.SetParam("MaxIter", (unsigned int)1000);
    solver_options.SetParam("PrintLevel", 0);

    hephaestus::InputParameters params;
    params.SetParam("UseGLVis", true);

    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("Order", 2);
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("DomainProperties", domain_properties);
    params.SetParam("FESpaces", fespaces);
    params.SetParam("GridFunctions", gridfunctions);
    params.SetParam("AuxKernels", auxkernels);
    params.SetParam("Postprocessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);
    params.SetParam("SolverOptions", solver_options);
    params.SetParam("Formulation", formulation);

    return params;
  }
};

TEST_F(TestComplexIrisWaveguide, CheckRun) {
  hephaestus::InputParameters params(test_params());
  hephaestus::SteadyExecutioner *executioner =
      new hephaestus::SteadyExecutioner(params);
  executioner->Init();
  executioner->Solve();
}
