#pragma once
#include "kernel_base.hpp"
#include "DomainLFH1DivIntegrator.hpp"

namespace hephaestus {

/*
(∇∙v, f)
*/
class DomainDivergenceLFKernel : public Kernel<mfem::ParLinearForm> {
public:
  DomainDivergenceLFKernel(const hephaestus::InputParameters &params);
  ~DomainDivergenceLFKernel();
  virtual void Init(hephaestus::GridFunctions &gridfunctions,
                    const hephaestus::FESpaces &fespaces,
                    hephaestus::BCMap &bc_map,
                    hephaestus::Coefficients &coefficients) override;
  virtual void Apply(mfem::ParLinearForm *lf) override;

  std::string coef_name;
  mfem::Coefficient *coef;
  
};

}; // namespace hephaestus
