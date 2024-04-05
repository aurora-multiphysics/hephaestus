#pragma once
#include "kernel_base.hpp"

namespace hephaestus
{

/*
(σ ∇ V, ∇ V')
*/
class DiffusionKernel : public Kernel<mfem::ParNonlinearForm>
{
public:
  DiffusionKernel(const hephaestus::InputParameters & params);

  ~DiffusionKernel() override = default;

  void Init(hephaestus::GridFunctions & gridfunctions,
            const hephaestus::FESpaces & fespaces,
            hephaestus::BCMap & bc_map,
            hephaestus::Coefficients & coefficients) override;
  void Apply(mfem::ParNonlinearForm * nlf) override;

  std::string _coef_name;
  mfem::Coefficient * _coef{nullptr};
};

} // namespace hephaestus
