#pragma once
#include "frequency_domain_formulation.hpp"

namespace hephaestus {
/*
Operator for solving:
∇×(α∇×u) + iωβu - ω²γu = g

via the weak form:
(α∇×u, ∇×u') + (iωβu, u') - (ω²ζu, u') - <(α∇×u)×n, u'> = (g, u')

where
u ∈ H(curl)
g ∈ H(div) is a divergence-free source field
ω is the angular frequency
α is the stiffness coefficient
ωβ is the loss coefficient
ω²γ is the mass coefficient

Dirichlet boundaries strongly constrain n×n×u
Integrated boundaries weakly constrain (α∇×u)×n
Robin boundaries weakly constrain (α∇×u)×n + γ(n×n×u) = F

Divergence cleaning (such as via Helmholtz projection)
should be performed on g before use in this operator.
*/
class ComplexMaxwellOperator : public EquationSystemOperator {
public:
  ComplexMaxwellOperator(mfem::ParMesh &pmesh, hephaestus::FESpaces &fespaces,
                         hephaestus::GridFunctions &gridfunctions,
                         hephaestus::BCMap &bc_map,
                         hephaestus::Coefficients &coefficients,
                         hephaestus::Sources &sources,
                         hephaestus::InputParameters &solver_options);

  ~ComplexMaxwellOperator(){};

  virtual void SetGridFunctions() override;
  virtual void Init(mfem::Vector &X) override;
  virtual void Solve(mfem::Vector &X) override;

  std::string h_curl_var_name, stiffness_coef_name, mass_coef_name,
      loss_coef_name;
  mfem::ParComplexGridFunction *u_;
  mfem::ParComplexLinearForm *b1_;
  mfem::ParSesquilinearForm *a1_;
  mfem::ComplexOperator::Convention conv_ = mfem::ComplexOperator::HERMITIAN;

  mfem::Coefficient *stiffCoef_; // Dia/Paramagnetic Material Coefficient
  mfem::Coefficient *massCoef_;  // -omega^2 epsilon
  mfem::Coefficient *lossCoef_;  // omega sigma

  mfem::VectorCoefficient *jrCoef_; // Volume Current Density Function
  mfem::VectorCoefficient *jiCoef_; // Volume Current Density Function

  mfem::Array<int> ess_bdr_tdofs_;
};

//
// Specifies output interfaces of a time-domain EM formulation.
class ComplexMaxwellFormulation
    : public hephaestus::FrequencyDomainFormulation {
  // std::vector<mfem::ParGridFunction *> local_trial_vars, local_test_vars;
protected:
  const std::string _alpha_coef_name;
  const std::string _beta_coef_name;
  const std::string _zeta_coef_name;
  const std::string _frequency_coef_name;
  const std::string _h_curl_var_name;
  const std::string _mass_coef_name;
  const std::string _loss_coef_name;

public:
  ComplexMaxwellFormulation(const std::string &frequency_coef_name,
                            const std::string &alpha_coef_name,
                            const std::string &beta_coef_name,
                            const std::string &zeta_coef_name,
                            const std::string &h_curl_var_name);

  virtual void ConstructOperator() override;

  virtual void RegisterGridFunctions() override;

  virtual void RegisterCoefficients() override;
};
} // namespace hephaestus
