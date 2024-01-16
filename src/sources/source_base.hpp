#pragma once
#include "../common/pfem_extras.hpp"
#include "coefficients.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "kernels.hpp"

namespace hephaestus
{

class Source : public hephaestus::Kernel<mfem::ParLinearForm>
{
public:
  Source() {}
  virtual void Init(hephaestus::GridFunctions & gridfunctions,
                    const hephaestus::FESpaces & fespaces,
                    hephaestus::BCMap & bc_map,
                    hephaestus::Coefficients & coefficients){};
  virtual void Apply(mfem::ParLinearForm * lf) override = 0;
  virtual void SubtractSource(mfem::ParGridFunction * gf) = 0;
};

}; // namespace hephaestus
