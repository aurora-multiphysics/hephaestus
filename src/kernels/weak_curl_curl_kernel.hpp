#pragma once
#include "kernel_base.hpp"

namespace hephaestus {

/*
(α∇×u_{n}, ∇×u')
*/
class WeakCurlCurlKernel : public Kernel<mfem::ParLinearForm> {
public:
  WeakCurlCurlKernel(const hephaestus::InputParameters &params);
  ~WeakCurlCurlKernel();
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::Coefficients &coefficients) override;
  virtual void Apply(mfem::ParLinearForm *lf) override;

  std::string coupled_gf_name;
  std::string coef_name;
  mfem::ParGridFunction *u_; //
  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::Coefficient *coef;
  mfem::ParBilinearForm *curlCurl;
};

}; // namespace hephaestus
