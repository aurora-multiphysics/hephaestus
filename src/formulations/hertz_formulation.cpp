#include "hertz_formulation.hpp"

namespace hephaestus {

HertzOperator::HertzOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : FrequencyDomainOperator(pmesh, order, fespaces, variables, bc_map,
                              domain_properties, sources, solver_options) {}

void HertzOperator::SetVariables() {
  state_var_names.push_back("electric_field_real");
  state_var_names.push_back("electric_field_imag");

  FrequencyDomainOperator::SetVariables();
};

void HertzOperator::Init(mfem::Vector &X) {
  FrequencyDomainOperator::Init(X);

  muInvCoef_ = _domain_properties.scalar_property_map["magnetic_reluctivity"];
  massCoef_ = _domain_properties.scalar_property_map["hertz_mass"];
  lossCoef_ = _domain_properties.scalar_property_map["hertz_loss"];
}

void HertzOperator::Solve(mfem::Vector &X) {
  mfem::OperatorHandle A1;
  mfem::Vector E, RHS;
  mfem::OperatorHandle PCOp;

  mfem::Vector zeroVec(3);
  zeroVec = 0.0;
  mfem::VectorConstantCoefficient zeroCoef(zeroVec);

  a1_ = new mfem::ParSesquilinearForm(e_->ParFESpace(), conv_);
  a1_->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*muInvCoef_), NULL);
  a1_->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*massCoef_), NULL);
  if (lossCoef_) {
    a1_->AddDomainIntegrator(NULL,
                             new mfem::VectorFEMassIntegrator(*lossCoef_));
  }

  // Volume Current Density
  mfem::Vector j(3);
  j = 0.0;
  jrCoef_ = new mfem::VectorConstantCoefficient(j);
  jiCoef_ = new mfem::VectorConstantCoefficient(j);

  jd_ = new mfem::ParComplexLinearForm(e_->ParFESpace(), conv_);
  jd_->AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*jrCoef_),
                           new mfem::VectorFEDomainLFIntegrator(*jiCoef_));

  _bc_map.applyEssentialBCs(std::string("electric_field"), ess_bdr_tdofs_, *e_,
                            pmesh_);
  _bc_map.applyIntegratedBCs(std::string("electric_field"), *jd_, pmesh_);
  _bc_map.applyIntegratedBCs(std::string("electric_field"), *a1_, pmesh_);

  a1_->Assemble();
  a1_->Finalize();

  jd_->Assemble();

  a1_->FormLinearSystem(ess_bdr_tdofs_, *e_, *jd_, A1, E, RHS);

  mfem::ComplexHypreParMatrix *A1Z = A1.As<mfem::ComplexHypreParMatrix>();
  mfem::HypreParMatrix *A1C = A1Z->GetSystemMatrix();
  mfem::SuperLURowLocMatrix A_SuperLU(*A1C);
  mfem::SuperLUSolver solver(MPI_COMM_WORLD);
  solver.SetOperator(A_SuperLU);
  solver.Mult(RHS, E);
  delete A1C;

  e_->Distribute(E);

  *_variables.Get(state_var_names.at(0)) = e_->real();
  *_variables.Get(state_var_names.at(1)) = e_->imag();
}

HertzFormulation::HertzFormulation() : SteadyFormulation() {
  frequency_coef_name = std::string("frequency");
  h_curl_var_name = std::string("electric_field");

  permittivity_coef_name = std::string("dielectric_permittivity");
  reluctivity_coef_name = std::string("magnetic_reluctivity");
  conductivity_coef_name = std::string("electrical_conductivity");
};

hephaestus::FrequencyDomainOperator *
HertzFormulation::CreateFrequencyDomainOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options) {
  fd_operator =
      new hephaestus::HertzOperator(pmesh, order, fespaces, variables, bc_map,
                                    domain_properties, sources, solver_options);
  return fd_operator;
};

void HertzFormulation::RegisterMissingVariables(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) {
  int myid;
  MPI_Comm_rank(pmesh.GetComm(), &myid);
  // Register default ParGridFunctions of state variables if not provided
  if (!variables.Has(h_curl_var_name)) {
    if (myid == 0) {
      std::cout
          << h_curl_var_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    fespaces.Register(
        "_HCurlFESpace",
        new mfem::common::ND_ParFESpace(&pmesh, 1, pmesh.Dimension()), true);
    variables.Register(h_curl_var_name + "_real",
                       new mfem::ParGridFunction(fespaces.Get("_HCurlFESpace")),
                       true);
    variables.Register(h_curl_var_name + "_imag",
                       new mfem::ParGridFunction(fespaces.Get("_HCurlFESpace")),
                       true);
  };
};

void HertzFormulation::RegisterCoefficients(
    hephaestus::DomainProperties &domain_properties) {

  freqCoef = dynamic_cast<mfem::ConstantCoefficient *>(
      domain_properties.scalar_property_map[frequency_coef_name]);
  if (freqCoef == NULL) {
    MFEM_ABORT("No frequency coefficient found. Frequency must be specified "
               "for frequency domain formulations.");
  }
  // define transformed
  domain_properties.scalar_property_map["_angular_frequency"] =
      new mfem::ConstantCoefficient(2.0 * M_PI * freqCoef->constant);
  domain_properties.scalar_property_map["_neg_angular_frequency"] =
      new mfem::ConstantCoefficient(-2.0 * M_PI * freqCoef->constant);
  domain_properties.scalar_property_map["_angular_frequency_sq"] =
      new mfem::ConstantCoefficient(pow(2.0 * M_PI * freqCoef->constant, 2));
  domain_properties.scalar_property_map["_neg_angular_frequency_sq"] =
      new mfem::ConstantCoefficient(-pow(2.0 * M_PI * freqCoef->constant, 2));

  if (domain_properties.scalar_property_map.count("magnetic_permeability") ==
      0) {
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
            std::string("magnetic_permeability")));
  }
  if (domain_properties.scalar_property_map.count(conductivity_coef_name) ==
      0) {
    domain_properties.scalar_property_map[conductivity_coef_name] =
        new mfem::PWCoefficient(
            domain_properties.getGlobalScalarProperty(conductivity_coef_name));
  }
  if (domain_properties.scalar_property_map.count(permittivity_coef_name) ==
      0) {
    domain_properties.scalar_property_map[permittivity_coef_name] =
        new mfem::PWCoefficient(
            domain_properties.getGlobalScalarProperty(permittivity_coef_name));
  }

  domain_properties.scalar_property_map["hertz_mass"] =
      new mfem::TransformedCoefficient(
          domain_properties.scalar_property_map["_neg_angular_frequency_sq"],
          domain_properties.scalar_property_map[permittivity_coef_name],
          prodFunc);

  domain_properties.scalar_property_map["hertz_loss"] =
      new mfem::TransformedCoefficient(
          domain_properties.scalar_property_map["_angular_frequency"],
          domain_properties.scalar_property_map[conductivity_coef_name],
          prodFunc);

  domain_properties.scalar_property_map[reluctivity_coef_name] =
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map["magnetic_permeability"],
          fracFunc);
}

} // namespace hephaestus
