#pragma once
#include "kernel_base.hpp"

namespace hephaestus {

/*
(σ ∇ V, u')
*/
class MixedVectorGradientKernel : public Kernel<mfem::ParMixedBilinearForm> {
public:
  MixedVectorGradientKernel(const hephaestus::InputParameters &params);
  virtual void Init(hephaestus::GridFunctions &gridfunctions,
                    const hephaestus::FESpaces &fespaces,
                    hephaestus::BCMap &bc_map,
                    hephaestus::Coefficients &coefficients) override;
  virtual void Apply(mfem::ParMixedBilinearForm *mblf) override;
  std::string coef_name;
  mfem::Coefficient *coef;
};

}; // namespace hephaestus
