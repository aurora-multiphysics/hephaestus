#pragma once
#include "../common/pfem_extras.hpp"
#include "auxkernels.hpp"
#include "equation_system.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

template <typename T>
std::vector<T *>
populateVectorFromNamedFieldsMap(mfem::NamedFieldsMap<T> nfmap,
                                 std::vector<std::string> keys) {
  std::vector<T *> result;
  for (auto &key : keys) {
    result.push_back(nfmap.Get(key));
  }
  return result;
};

// Specifies output interfaces of a time-domain EM formulation.
class TransientFormulation : public mfem::TimeDependentOperator {

public:
  TransientFormulation(){};

  ~TransientFormulation(){};

  virtual void Init(mfem::Vector &X){};

  virtual void RegisterVariables() = 0;
  virtual void RegisterAuxKernels(hephaestus::AuxKernels &auxkernels){};
  std::string GetTimeDerivativeName(std::string name) {
    return std::string("d") + name + std::string("_dt");
  }
  std::vector<std::string>
  GetTimeDerivativeNames(std::vector<std::string> gridfunction_names) {
    std::vector<std::string> time_derivative_names;
    for (auto &gridfunction_name : gridfunction_names) {
      time_derivative_names.push_back(GetTimeDerivativeName(gridfunction_name));
    }
    return time_derivative_names;
  }
  std::vector<mfem::ParGridFunction *> registerTimeDerivatives(
      std::vector<std::string> gridfunction_names,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &gridfunctions) {
    std::vector<mfem::ParGridFunction *> time_derivatives;

    for (auto &gridfunction_name : gridfunction_names) {
      gridfunctions.Register(
          GetTimeDerivativeName(gridfunction_name),
          new mfem::ParGridFunction(
              gridfunctions.Get(gridfunction_name)->ParFESpace()),
          true);
      time_derivatives.push_back(
          gridfunctions.Get(GetTimeDerivativeName(gridfunction_name)));
    }
    return time_derivatives;
  }
  mfem::Array<int> true_offsets;
  // Vector of names of state variables used in formulation, ordered by
  // appearance in block vector during solve.
  std::vector<std::string> state_var_names;
  // Vector of names of recognised auxiliary variables that can be calculated
  // from formulation,
  std::vector<std::string> aux_var_names;
  // Vector of names of active auxiliary variables that are being calculated
  // in formulation,
  std::vector<std::string> active_aux_var_names;

  std::vector<mfem::ParGridFunction *> local_trial_vars, local_test_vars;

  mfem::ConstantCoefficient oneCoef =
      mfem::ConstantCoefficient(1.0); // Auxiliary coefficient
};
} // namespace hephaestus
