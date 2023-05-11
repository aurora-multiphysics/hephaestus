#include "complex_maxwell_form.hpp"

namespace hephaestus {

ComplexMaxwellOperator::ComplexMaxwellOperator(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : FrequencyDomainOperator(pmesh, fespaces, variables, bc_map,
                              domain_properties, sources, solver_options),
      h_curl_var_name(solver_options.GetParam<std::string>("HCurlVarName")),
      stiffness_coef_name(
          solver_options.GetParam<std::string>("StiffnessCoefName")),
      mass_coef_name(solver_options.GetParam<std::string>("MassCoefName")),
      loss_coef_name(solver_options.GetParam<std::string>("LossCoefName")) {}

void ComplexMaxwellOperator::SetVariables() {
  state_var_names.push_back(h_curl_var_name + "_real");
  state_var_names.push_back(h_curl_var_name + "_imag");

  FrequencyDomainOperator::SetVariables();

  u_ = new mfem::ParComplexGridFunction(local_test_vars.at(0)->ParFESpace());
  *u_ = std::complex(0.0, 0.0);
};

void ComplexMaxwellOperator::Init(mfem::Vector &X) {
  FrequencyDomainOperator::Init(X);

  stiffCoef_ = _domain_properties.scalar_property_map[stiffness_coef_name];
  massCoef_ = _domain_properties.scalar_property_map[mass_coef_name];
  lossCoef_ = _domain_properties.scalar_property_map[loss_coef_name];
}

void ComplexMaxwellOperator::Solve(mfem::Vector &X) {
  mfem::OperatorHandle A1;
  mfem::Vector U, RHS;
  mfem::OperatorHandle PCOp;

  mfem::Vector zeroVec(3);
  zeroVec = 0.0;
  mfem::VectorConstantCoefficient zeroCoef(zeroVec);

  a1_ = new mfem::ParSesquilinearForm(u_->ParFESpace(), conv_);
  a1_->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*stiffCoef_), NULL);
  if (massCoef_) {
    a1_->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*massCoef_),
                             NULL);
  }
  if (lossCoef_) {
    a1_->AddDomainIntegrator(NULL,
                             new mfem::VectorFEMassIntegrator(*lossCoef_));
  }

  // Volume Current Density
  mfem::Vector j(3);
  j = 0.0;
  jrCoef_ = new mfem::VectorConstantCoefficient(j);
  jiCoef_ = new mfem::VectorConstantCoefficient(j);

  mfem::ParLinearForm *b1_real_ = new mfem::ParLinearForm(u_->ParFESpace());
  mfem::ParLinearForm *b1_imag_ = new mfem::ParLinearForm(u_->ParFESpace());
  *b1_real_ = 0.0;
  *b1_imag_ = 0.0;

  _sources.Apply(b1_real_);

  b1_ = new mfem::ParComplexLinearForm(u_->ParFESpace(), conv_);

  _bc_map.applyEssentialBCs(std::string(h_curl_var_name), ess_bdr_tdofs_, *u_,
                            pmesh_);
  _bc_map.applyIntegratedBCs(std::string(h_curl_var_name), *b1_, pmesh_);
  _bc_map.applyIntegratedBCs(std::string(h_curl_var_name), *a1_, pmesh_);

  a1_->Assemble();
  a1_->Finalize();

  b1_->Assemble();
  b1_->real() += *b1_real_;
  b1_->imag() += *b1_imag_;

  a1_->FormLinearSystem(ess_bdr_tdofs_, *u_, *b1_, A1, U, RHS);

  mfem::ComplexHypreParMatrix *A1Z = A1.As<mfem::ComplexHypreParMatrix>();
  mfem::HypreParMatrix *A1C = A1Z->GetSystemMatrix();
  mfem::SuperLURowLocMatrix A_SuperLU(*A1C);
  mfem::SuperLUSolver solver(MPI_COMM_WORLD);
  solver.SetOperator(A_SuperLU);
  solver.Mult(RHS, U);
  delete A1C;

  a1_->RecoverFEMSolution(U, *b1_, *u_);

  *_variables.Get(state_var_names.at(0)) = u_->real();
  *_variables.Get(state_var_names.at(1)) = u_->imag();
}

ComplexMaxwellFormulation::ComplexMaxwellFormulation() : SteadyFormulation() {
  frequency_coef_name = std::string("frequency");
  h_curl_var_name = std::string("h_curl_var");
  mass_coef_name = std::string("maxwell_mass");
  loss_coef_name = std::string("maxwell_loss");

  zeta_coef_name = std::string("zeta");
  alpha_coef_name = std::string("alpha");
  beta_coef_name = std::string("beta");
};

hephaestus::FrequencyDomainOperator *
ComplexMaxwellFormulation::CreateFrequencyDomainOperator(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options) {

  solver_options.SetParam("HCurlVarName", h_curl_var_name);
  solver_options.SetParam("StiffnessCoefName", alpha_coef_name);
  solver_options.SetParam("MassCoefName", mass_coef_name);
  solver_options.SetParam("LossCoefName", loss_coef_name);

  fd_operator = new hephaestus::ComplexMaxwellOperator(
      pmesh, fespaces, variables, bc_map, domain_properties, sources,
      solver_options);
  return fd_operator;
};

void ComplexMaxwellFormulation::RegisterMissingVariables(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) {
  int myid;
  MPI_Comm_rank(pmesh.GetComm(), &myid);
  // Register default ParGridFunctions of state variables if not provided
  if (!variables.Has(h_curl_var_name + "_real")) {
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

void ComplexMaxwellFormulation::RegisterCoefficients(
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
  if (domain_properties.scalar_property_map.count(beta_coef_name) == 0) {
    domain_properties.scalar_property_map[beta_coef_name] =
        new mfem::PWCoefficient(
            domain_properties.getGlobalScalarProperty(beta_coef_name));
  }
  if (domain_properties.scalar_property_map.count(zeta_coef_name) == 0) {
    domain_properties.scalar_property_map[zeta_coef_name] =
        new mfem::PWCoefficient(
            domain_properties.getGlobalScalarProperty(zeta_coef_name));
  }

  domain_properties.scalar_property_map[mass_coef_name] =
      new mfem::TransformedCoefficient(
          domain_properties.scalar_property_map["_neg_angular_frequency_sq"],
          domain_properties.scalar_property_map[zeta_coef_name], prodFunc);

  domain_properties.scalar_property_map[loss_coef_name] =
      new mfem::TransformedCoefficient(
          domain_properties.scalar_property_map["_angular_frequency"],
          domain_properties.scalar_property_map[beta_coef_name], prodFunc);

  domain_properties.scalar_property_map[alpha_coef_name] =
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map["magnetic_permeability"],
          fracFunc);
}

} // namespace hephaestus