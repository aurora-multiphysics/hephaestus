#pragma once
#include "kernel_base.hpp"
#include <set>
#include <memory>
#include <mfem.hpp>

namespace hephaestus
{
/// Lightweight store for kernels of a particular type. Ensures that kernels are unique and non-null.
template <typename TKernel>
class KernelStore
{
public:
  using SetType = std::set<std::shared_ptr<TKernel>>;
  using iterator = typename SetType::iterator;

  KernelStore() = default;
  ~KernelStore() = default;

  /// Add a new kernel to the store.
  inline void Register(std::shared_ptr<TKernel> kernel)
  {
    CheckKernelIsRegistrable(kernel);

    _kernels.insert(std::move(kernel));
  }

  /// Remove a kernel from the store.
  inline void Deregister(const std::shared_ptr<TKernel> & kernel) { _kernels.erase(kernel); }

  /// Remove all kernels from the store.
  inline void DeregisterAll() { _kernels.clear(); }

  // NOLINTNEXTLINE(readability-identifier-naming)
  inline iterator begin() { return _kernels.begin(); }

  // NOLINTNEXTLINE(readability-identifier-naming)
  inline iterator end() { return _kernels.end(); }

protected:
  /// Ensure that the kernel is non-null and is not already registered.
  void CheckKernelIsRegistrable(std::shared_ptr<TKernel> & kernel) const
  {
    if (!kernel)
    {
      MFEM_ABORT("Cannot register NULL kernel.");
    }

    if (_kernels.count(kernel) > 0)
    {
      MFEM_ABORT("The kernel is already registered. Cannot register twice.");
    }
  }

private:
  std::set<std::shared_ptr<TKernel>> _kernels;
};

// Typedefs.
using ParBilinearFormKernelStore = KernelStore<ParBilinearFormKernel>;
using ParLinearFormKernelStore = KernelStore<ParLinearFormKernel>;
using ParNonlinearFormKernelStore = KernelStore<ParNonlinearFormKernel>;
using ParMixedBilinearFormKernelStore = KernelStore<ParMixedBilinearFormKernel>;
}
