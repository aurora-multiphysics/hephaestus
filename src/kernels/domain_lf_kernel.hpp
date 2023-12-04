#pragma once
#include "kernel_base.hpp"

namespace hephaestus {

/*
(σ ∇ V, ∇ V')
*/
class DomainLFKernel : public Kernel<mfem::ParLinearForm> {
public:
  DomainLFKernel(const hephaestus::InputParameters &params);
  virtual void Init(hephaestus::GridFunctions &gridfunctions,
                    const hephaestus::FESpaces &fespaces,
                    hephaestus::BCMap &bc_map,
                    hephaestus::Coefficients &coefficients) override;
  virtual void Apply(mfem::ParLinearForm *lf) override;
  std::string coef_name;
  mfem::Coefficient *coef;
};

}; // namespace hephaestus
