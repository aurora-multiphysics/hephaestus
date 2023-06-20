#pragma once
#include "kernel_base.hpp"

namespace hephaestus {

/*
(σ ∇ V, u')
*/
class MixedVectorGradientKernel : public Kernel<mfem::ParMixedBilinearForm> {
public:
  MixedVectorGradientKernel(const hephaestus::InputParameters &params);
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map,
       hephaestus::DomainProperties &domain_properties) override;
  virtual void Apply(mfem::ParMixedBilinearForm *mblf) override;
  std::string coef_name;
  mfem::Coefficient *coef;
};

}; // namespace hephaestus
