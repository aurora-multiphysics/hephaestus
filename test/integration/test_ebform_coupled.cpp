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

class JouleHeatingCoefficient : public hephaestus::CoupledCoefficient {
private:
  mfem::Coefficient *sigma;
  std::string conductivity_coef_name;

public:
  JouleHeatingCoefficient(const hephaestus::InputParameters &params)
      : hephaestus::CoupledCoefficient(params),
        conductivity_coef_name(
            params.GetParam<std::string>("ConductivityCoefName")){};

  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::DomainProperties &domain_properties) {
    hephaestus::CoupledCoefficient::Init(variables, domain_properties);
    sigma = domain_properties.scalar_property_map[conductivity_coef_name];
  }
  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip) {
    mfem::Vector E;
    double thisSigma;
    gf->GetVectorValue(T, ip, E);
    thisSigma = sigma->Eval(T, ip);
    return thisSigma * (E * E);
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

    hephaestus::FESpaces fespaces;
    hephaestus::GridFunctions gridfunctions;

    hephaestus::InputParameters l2fespaceparams;
    l2fespaceparams.SetParam("FESpaceName", std::string("L2"));
    l2fespaceparams.SetParam("FESpaceType", std::string("L2"));
    l2fespaceparams.SetParam("order", 1);
    l2fespaceparams.SetParam("components", 3);
    fespaces.StoreInput(l2fespaceparams);

    hephaestus::InputParameters jouleheatinggfparams;
    jouleheatinggfparams.SetParam("VariableName",
                                  std::string("joule_heating_load"));
    jouleheatinggfparams.SetParam("FESpaceName", std::string("L2"));
    gridfunctions.StoreInput(jouleheatinggfparams);

    hephaestus::InputParameters h1fespaceparams;
    h1fespaceparams.SetParam("FESpaceName", std::string("H1"));
    h1fespaceparams.SetParam("FESpaceType", std::string("H1"));
    h1fespaceparams.SetParam("order", 1);
    h1fespaceparams.SetParam("components", 3);
    fespaces.StoreInput(h1fespaceparams);

    hephaestus::InputParameters hcurlsourcefespaceparams;
    hcurlsourcefespaceparams.SetParam("FESpaceName",
                                      std::string("HCurlSource"));
    hcurlsourcefespaceparams.SetParam("FESpaceType", std::string("ND"));
    hcurlsourcefespaceparams.SetParam("order", 2);
    hcurlsourcefespaceparams.SetParam("components", 3);
    hephaestus::InputParameters h1sourcefespaceparams;
    h1sourcefespaceparams.SetParam("FESpaceName", std::string("H1Source"));
    h1sourcefespaceparams.SetParam("FESpaceType", std::string("H1"));
    h1sourcefespaceparams.SetParam("order", 2);
    h1sourcefespaceparams.SetParam("components", 3);
    fespaces.StoreInput(hcurlsourcefespaceparams);
    fespaces.StoreInput(h1sourcefespaceparams);

    hephaestus::InputParameters temperaturegfparams;
    temperaturegfparams.SetParam("VariableName", std::string("temperature"));
    temperaturegfparams.SetParam("FESpaceName", std::string("H1"));
    gridfunctions.StoreInput(temperaturegfparams);

    // materialCopper instances and get property coefs? init can be for all...
    // CoupledCoefficients must also be added to Auxkernels
    hephaestus::InputParameters copper_conductivity_params;
    copper_conductivity_params.SetParam("CoupledVariableName",
                                        std::string("temperature"));
    hephaestus::CoupledCoefficient *wireConductivity =
        new CopperConductivityCoefficient(copper_conductivity_params);

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

    hephaestus::InputParameters joule_heating_params;
    joule_heating_params.SetParam("CoupledVariableName",
                                  std::string("electric_field"));
    joule_heating_params.SetParam("ConductivityCoefName",
                                  std::string("electrical_conductivity"));
    JouleHeatingCoefficient *jouleHeating =
        new JouleHeatingCoefficient(joule_heating_params);
    domain_properties.scalar_property_map["JouleHeating"] = jouleHeating;

    hephaestus::AuxKernels auxkernels;
    auxkernels.Register("CoupledCoefficient", wireConductivity, false);
    auxkernels.Register("JouleHeating", jouleHeating, false);

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

    mfem::Mesh mesh(
        (std::string(DATA_DIR) + std::string("./cylinder-hex-q2.gen")).c_str(),
        1, 1);

    std::map<std::string, mfem::DataCollection *> data_collections;
    data_collections["VisItDataCollection"] =
        new mfem::VisItDataCollection("EBFormVisIt");
    data_collections["ParaViewDataCollection"] =
        new mfem::ParaViewDataCollection("EBFormParaView");
    hephaestus::Outputs outputs(data_collections);

    hephaestus::InputParameters jouleheatingauxparams;
    jouleheatingauxparams.SetParam("VariableName",
                                   std::string("joule_heating_load"));
    jouleheatingauxparams.SetParam("CoefficientName",
                                   std::string("JouleHeating"));
    auxkernels.Register(
        "CoefficientAuxKernel",
        new hephaestus::CoefficientAuxKernel(jouleheatingauxparams), true);

    hephaestus::Postprocessors postprocessors;
    hephaestus::Sources sources;
    hephaestus::InputParameters scalar_potential_source_params;
    scalar_potential_source_params.SetParam("SourceName",
                                            std::string("source"));
    scalar_potential_source_params.SetParam("PotentialName",
                                            std::string("electric_potential"));
    scalar_potential_source_params.SetParam("HCurlFESpaceName",
                                            std::string("HCurlSource"));
    scalar_potential_source_params.SetParam("H1FESpaceName",
                                            std::string("H1Source"));
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
    params.SetParam("Order", 2);
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

TEST_F(TestEBFormCoupled, CheckRun) {
  hephaestus::InputParameters params(test_params());
  hephaestus::TransientExecutioner *executioner =
      new hephaestus::TransientExecutioner(params);
  executioner->Init();
  executioner->Solve();
}
