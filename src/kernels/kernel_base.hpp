#pragma once
#include "../common/pfem_extras.hpp"
#include "coefficients.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"

namespace hephaestus {

template <typename T> class Kernel {
public:
  Kernel() {}
  Kernel(const hephaestus::InputParameters &params){};
  virtual void
  Init(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
       hephaestus::BCMap &bc_map, hephaestus::Coefficients &coefficients){};

  virtual void Apply(T *form) = 0;
};

}; // namespace hephaestus
