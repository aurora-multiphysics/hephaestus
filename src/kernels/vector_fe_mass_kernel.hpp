#pragma once
#include "kernel_base.hpp"

namespace hephaestus
{

/*
(βu, u')
*/
class VectorFEMassKernel : public Kernel<mfem::ParNonlinearForm>
{
public:
  VectorFEMassKernel(const hephaestus::InputParameters & params);

  ~VectorFEMassKernel() override = default;

  void Init(hephaestus::GridFunctions & gridfunctions,
            const hephaestus::FESpaces & fespaces,
            hephaestus::BCMap & bc_map,
            hephaestus::Coefficients & coefficients) override;
  void Apply(mfem::ParNonlinearForm * nlf) override;
  std::string _coef_name;
  mfem::Coefficient * _coef{nullptr};
};

} // namespace hephaestus
