#pragma once
#include "kernel_base.hpp"

namespace hephaestus {

/*
(αv_{n}, ∇×u')
*/
class WeakCurlKernel : public Kernel<mfem::ParLinearForm> {
public:
  WeakCurlKernel(const hephaestus::InputParameters &params);
  ~WeakCurlKernel();
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::Coefficients &domain_properties) override;
  virtual void Apply(mfem::ParLinearForm *lf) override;

  std::string hcurl_gf_name, hdiv_gf_name;
  std::string coef_name;
  mfem::ParGridFunction *u_, *v_; //
  mfem::Coefficient *coef;
  mfem::ParMixedBilinearForm *weakCurl;
};

}; // namespace hephaestus
