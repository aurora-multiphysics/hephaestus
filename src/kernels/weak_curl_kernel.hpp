#pragma once
#include "kernel_base.hpp"

namespace hephaestus {

/*
(αv_{n}, ∇×u')
*/
class WeakCurlKernel : public Kernel<mfem::ParLinearForm> {
public:
  WeakCurlKernel(const hephaestus::InputParameters &params);

  ~WeakCurlKernel() override{};

  virtual void Init(hephaestus::GridFunctions &gridfunctions,
                    const hephaestus::FESpaces &fespaces,
                    hephaestus::BCMap &bc_map,
                    hephaestus::Coefficients &coefficients) override;
  virtual void Apply(mfem::ParLinearForm *lf) override;

  std::string hcurl_gf_name, hdiv_gf_name;
  std::string coef_name;
  mfem::ParGridFunction *u_, *v_; //
  mfem::Coefficient *coef;

  std::unique_ptr<mfem::ParMixedBilinearForm> weakCurl;
};

}; // namespace hephaestus
