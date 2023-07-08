#pragma once
#include "kernel_base.hpp"

namespace hephaestus {

/*
(σ u, ∇ V')
*/
class VectorFEWeakDivergenceKernel : public Kernel<mfem::ParMixedBilinearForm> {
public:
  VectorFEWeakDivergenceKernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::Coefficients &coefficients) override;
  virtual void Apply(mfem::ParMixedBilinearForm *mblf) override;
  std::string coef_name;
  mfem::Coefficient *coef;
};

}; // namespace hephaestus
