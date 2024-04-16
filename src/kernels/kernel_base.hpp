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
    _gridfunctions = &gridfunctions;
    _fespaces = &fespaces;
    _bc_map = &bc_map;
    _coefficients = &coefficients;
  }

  virtual void Apply(T * form) = 0;

  // Update method. All kernels should override this.
  virtual void Update()
  {
    // Call correct Init method again for kernel.
    Init(*const_cast<GridFunctions *>(_gridfunctions),
         *const_cast<FESpaces *>(_fespaces),
         *const_cast<BCMap *>(_bc_map),
         *const_cast<Coefficients *>(_coefficients));
  }

private: // NB: - temporary.
  const hephaestus::GridFunctions * _gridfunctions{nullptr};
  const hephaestus::FESpaces * _fespaces{nullptr};
  const hephaestus::BCMap * _bc_map{nullptr};
  const hephaestus::Coefficients * _coefficients{nullptr};
};

} // namespace hephaestus
