#pragma once
#include "steady_formulation.hpp"

namespace hephaestus {
// Specifies output interfaces of a time-domain EM formulation.
//   Curl mu^{-1} Curl A  + i omega sigma A = J
class ComplexAFormOperator : public FrequencyDomainOperator {
public:
  ComplexAFormOperator(
      mfem::ParMesh &pmesh, int order,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options);

  ~ComplexAFormOperator(){};

  virtual void SetVariables() override;
  virtual void Init(mfem::Vector &X) override;
  virtual void Solve(mfem::Vector &X) override;

  mfem::ParComplexGridFunction *a_;
  mfem::ParLinearForm *j_real_;
  mfem::ParLinearForm *j_imag_;
  mfem::ParComplexLinearForm *jd_;
  mfem::ParSesquilinearForm *a1_;
  mfem::ComplexOperator::Convention conv_ = mfem::ComplexOperator::HERMITIAN;

  mfem::Coefficient *muInvCoef_; // Dia/Paramagnetic Material Coefficient
  mfem::Coefficient *lossCoef_;  // omega sigma

  mfem::VectorCoefficient *jrCoef_; // Volume Current Density Function
  mfem::VectorCoefficient *jiCoef_; // Volume Current Density Function

  mfem::Array<int> ess_bdr_tdofs_;
};

//
// Specifies output interfaces of a time-domain EM formulation.
class ComplexAFormulation : public SteadyFormulation {
  // std::vector<mfem::ParGridFunction *> local_trial_vars, local_test_vars;
  std::string frequency_coef_name, permittivity_coef_name,
      reluctivity_coef_name, conductivity_coef_name, h_curl_var_name;

public:
  ComplexAFormulation();

  // virtual hephaestus::TimeDependentEquationSystem *CreateEquationSystem();

  virtual hephaestus::FrequencyDomainOperator *CreateFrequencyDomainOperator(
      mfem::ParMesh &pmesh, int order,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options) override;

  virtual void RegisterMissingVariables(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) override;

  virtual void
  RegisterAuxKernels(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                     hephaestus::AuxKernels &auxkernels) override;

  virtual void RegisterCoefficients(
      hephaestus::DomainProperties &domain_properties) override;
};
} // namespace hephaestus
