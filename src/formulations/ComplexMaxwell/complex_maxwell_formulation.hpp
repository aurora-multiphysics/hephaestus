#pragma once
#include "frequency_domain_em_formulation.hpp"

namespace hephaestus {
/*
Forumulation for solving:
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
class ComplexMaxwellFormulation
    : public hephaestus::FrequencyDomainEMFormulation {
protected:
  const std::string _alpha_coef_name;
  const std::string _beta_coef_name;
  const std::string _zeta_coef_name;
  const std::string _frequency_coef_name;
  const std::string _h_curl_var_complex_name;
  const std::string _h_curl_var_real_name;
  const std::string _h_curl_var_imag_name;
  const std::string _mass_coef_name;
  const std::string _loss_coef_name;

public:
  ComplexMaxwellFormulation(const std::string &frequency_coef_name,
                            const std::string &alpha_coef_name,
                            const std::string &beta_coef_name,
                            const std::string &zeta_coef_name,
                            const std::string &h_curl_var_complex_name,
                            const std::string &h_curl_var_real_name,
                            const std::string &h_curl_var_imag_name);

  virtual void ConstructJacobianSolver() override;

  virtual void ConstructOperator() override;

  virtual void RegisterGridFunctions() override;

  virtual void RegisterCoefficients() override;
};

class ComplexMaxwellOperator : public ProblemOperator {
public:
  ComplexMaxwellOperator(hephaestus::Problem &problem,
                         const std::string &h_curl_var_complex_name,
                         const std::string &h_curl_var_real_name,
                         const std::string &h_curl_var_imag_name,
                         const std::string &stiffness_coef_name,
                         const std::string &mass_coef_name,
                         const std::string &loss_coef_name);

  ~ComplexMaxwellOperator(){};

  virtual void SetGridFunctions() override;
  virtual void Init(mfem::Vector &X) override;
  virtual void Solve(mfem::Vector &X) override;

  std::string _h_curl_var_complex_name, _h_curl_var_real_name,
      _h_curl_var_imag_name, _stiffness_coef_name, _mass_coef_name,
      _loss_coef_name;
  mfem::ParComplexGridFunction *u_;
  mfem::ComplexOperator::Convention conv_ = mfem::ComplexOperator::HERMITIAN;

  mfem::Coefficient *stiffCoef_; // Dia/Paramagnetic Material Coefficient
  mfem::Coefficient *massCoef_;  // -omega^2 epsilon
  mfem::Coefficient *lossCoef_;  // omega sigma

  mfem::Array<int> ess_bdr_tdofs_;
};

} // namespace hephaestus
