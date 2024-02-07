#pragma once
#include "../common/pfem_extras.hpp"
#include "coefficients.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"

namespace hephaestus
{

template <typename T>
class Kernel
{
public:
  Kernel() = default;
  virtual ~Kernel() = default;

  Kernel(const hephaestus::InputParameters & params) {}
  virtual void Init(hephaestus::GridFunctions & gridfunctions,
                    const hephaestus::FESpaces & fespaces,
                    hephaestus::BCMap & bc_map,
                    hephaestus::Coefficients & coefficients)
  {
  }

  virtual void Apply(T * form) = 0;
};

// Typedefs.
using ParBilinearFormKernel = hephaestus::Kernel<mfem::ParBilinearForm>;
using ParLinearFormKernel = hephaestus::Kernel<mfem::ParLinearForm>;
using ParNonlinearFormKernel = hephaestus::Kernel<mfem::ParNonlinearForm>;
using ParMixedBilinearFormKernel = hephaestus::Kernel<mfem::ParMixedBilinearForm>;

} // namespace hephaestus
