#pragma once
#include "kernel_base.hpp"

namespace hephaestus {

/*
(α∇×u, ∇×u')
*/
class CurlCurlKernel : public Kernel<mfem::ParBilinearForm> {
public:
  CurlCurlKernel(const hephaestus::InputParameters &params);

  ~CurlCurlKernel() override{};

  virtual void Init(hephaestus::GridFunctions &gridfunctions,
                    const hephaestus::FESpaces &fespaces,
                    hephaestus::BCMap &bc_map,
                    hephaestus::Coefficients &coefficients) override;
  virtual void Apply(mfem::ParBilinearForm *blf) override;
  std::string coef_name;
  mfem::Coefficient *coef;
};

}; // namespace hephaestus
