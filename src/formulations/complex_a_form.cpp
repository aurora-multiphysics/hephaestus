#include "complex_a_form.hpp"

namespace hephaestus {

ComplexAFormOperator::ComplexAFormOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : FrequencyDomainOperator(pmesh, order, fespaces, variables, bc_map,
                              domain_properties, sources, solver_options) {}

void ComplexAFormOperator::SetVariables() {
  state_var_names.push_back("magnetic_vector_potential_real");
  state_var_names.push_back("magnetic_vector_potential_imag");

  FrequencyDomainOperator::SetVariables();

  a_ = new mfem::ParComplexGridFunction(local_test_vars.at(0)->ParFESpace());
  *a_ = std::complex(0.0, 0.0);
};

void ComplexAFormOperator::Init(mfem::Vector &X) {
  FrequencyDomainOperator::Init(X);

  muInvCoef_ = _domain_properties.scalar_property_map["magnetic_reluctivity"];
  lossCoef_ = _domain_properties.scalar_property_map["hertz_loss"];
}

void ComplexAFormOperator::Solve(mfem::Vector &X) {
  mfem::OperatorHandle A1;
  mfem::Vector A, RHS;
  mfem::OperatorHandle PCOp;

  mfem::Vector zeroVec(3);
  zeroVec = 0.0;
  mfem::VectorConstantCoefficient zeroCoef(zeroVec);
  a1_ = new mfem::ParSesquilinearForm(a_->ParFESpace(), conv_);
  a1_->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*muInvCoef_), NULL);
  if (lossCoef_) {
    a1_->AddDomainIntegrator(NULL,
                             new mfem::VectorFEMassIntegrator(*lossCoef_));
  }

  // Volume Current Density
  mfem::Vector j(3);
  j = 0.0;
  jrCoef_ = new mfem::VectorConstantCoefficient(j);
  jiCoef_ = new mfem::VectorConstantCoefficient(j);
  j_real_ = new mfem::ParLinearForm(a_->ParFESpace());
  j_imag_ = new mfem::ParLinearForm(a_->ParFESpace());

  *j_real_ = 0.0;
  *j_imag_ = 0.0;

  _sources.Apply(j_real_);
  // _sources.Apply(j_imag_);
  jd_ =
      new mfem::ParComplexLinearForm(a_->ParFESpace(), j_real_, j_imag_, conv_);

  _bc_map.applyEssentialBCs(std::string("magnetic_vector_potential"),
                            ess_bdr_tdofs_, *a_, pmesh_);
  _bc_map.applyIntegratedBCs(std::string("magnetic_vector_potential"), *jd_,
                             pmesh_);
  _bc_map.applyIntegratedBCs(std::string("magnetic_vector_potential"), *a1_,
                             pmesh_);

  a1_->Assemble();
  a1_->Finalize();

  // jd_->Assemble();
  jd_->real() = *j_real_;
  jd_->imag() = *j_imag_;

  a1_->FormLinearSystem(ess_bdr_tdofs_, *a_, *jd_, A1, A, RHS);

  mfem::ComplexHypreParMatrix *A1Z = A1.As<mfem::ComplexHypreParMatrix>();
  mfem::HypreParMatrix *A1C = A1Z->GetSystemMatrix();
  mfem::SuperLURowLocMatrix A_SuperLU(*A1C);
  mfem::SuperLUSolver solver(MPI_COMM_WORLD);
  solver.SetOperator(A_SuperLU);
  solver.Mult(RHS, A);
  delete A1C;

  a1_->RecoverFEMSolution(A, *jd_, *a_);

  *_variables.Get(state_var_names.at(0)) = a_->real();
  *_variables.Get(state_var_names.at(1)) = a_->imag();
}

ComplexAFormulation::ComplexAFormulation() : SteadyFormulation() {
  frequency_coef_name = std::string("frequency");
  h_curl_var_name = std::string("magnetic_vector_potential");

  permittivity_coef_name = std::string("dielectric_permittivity");
  reluctivity_coef_name = std::string("magnetic_reluctivity");
  conductivity_coef_name = std::string("electrical_conductivity");
};

hephaestus::FrequencyDomainOperator *
ComplexAFormulation::CreateFrequencyDomainOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options) {
  fd_operator = new hephaestus::ComplexAFormOperator(
      pmesh, order, fespaces, variables, bc_map, domain_properties, sources,
      solver_options);
  return fd_operator;
};

void ComplexAFormulation::RegisterMissingVariables(
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

void ComplexAFormulation::RegisterAuxKernels(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::AuxKernels &auxkernels) {
  std::vector<std::string> aux_var_names;
  std::string b_field_name = "magnetic_flux_density";
  if (variables.Get(b_field_name + "_real") != NULL) {
    // if (myid_ == 0) {
    std::cout << b_field_name + "_real"
              << " found in variables: building auxvar " << std::endl;
    // }
    hephaestus::InputParameters b_field_aux_params;
    b_field_aux_params.SetParam("VariableName", h_curl_var_name + "_real");
    b_field_aux_params.SetParam("CurlVariableName", b_field_name + "_real");
    auxkernels.Register("_magnetic_flux_density_re_aux",
                        new hephaestus::CurlAuxKernel(b_field_aux_params),
                        true);

    b_field_aux_params.SetParam("VariableName", h_curl_var_name + "_imag");
    b_field_aux_params.SetParam("CurlVariableName", b_field_name + "_imag");
    auxkernels.Register("_magnetic_flux_density_im_aux",
                        new hephaestus::CurlAuxKernel(b_field_aux_params),
                        true);
  }
}

void ComplexAFormulation::RegisterCoefficients(
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
