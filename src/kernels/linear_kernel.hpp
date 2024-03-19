#pragma once
#include "kernel_base.hpp"

namespace hephaestus
{

class LinearKernel : public Kernel<mfem::ParLinearForm>
{
public:
  LinearKernel(const hephaestus::InputParameters & params);

  ~LinearKernel() override = default;

  void Init(hephaestus::GridFunctions & gridfunctions,
            const hephaestus::FESpaces & fespaces,
            hephaestus::BCMap & bc_map,
            hephaestus::Coefficients & coefficients) override;
  void Apply(mfem::ParLinearForm * blf) override;

  std::string _coef_name;
  mfem::Coefficient * _coef{nullptr};
};

} // namespace hephaestus
