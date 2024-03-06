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

// Aliases.
using ParBilinearFormKernel = Kernel<mfem::ParBilinearForm>;
using ParLinearFormKernel = Kernel<mfem::ParLinearForm>;
using ParNonlinearFormKernel = Kernel<mfem::ParNonlinearForm>;
using ParMixedBilinearFormKernel = Kernel<mfem::ParMixedBilinearForm>;

template <typename T>
using KernelStore = std::vector<std::shared_ptr<T>>;

using ParBilinearFormKernelStore = KernelStore<ParBilinearFormKernel>;
using ParLinearFormKernelStore = KernelStore<ParLinearFormKernel>;
using ParNonlinearFormKernelStore = KernelStore<ParNonlinearFormKernel>;
using ParMixedBilinearFormKernelStore = KernelStore<ParMixedBilinearFormKernel>;

} // namespace hephaestus
