#pragma once
#include "../common/pfem_extras.hpp"
#include "factory.hpp"
#include "inputs.hpp"
#include "mfem.hpp"

namespace hephaestus {

class VariableMap : public mfem::NamedFieldsMap<mfem::ParGridFunction> {};

class Variables {
public:
  Variables() {}

  Variables(const hephaestus::InputParameters);

  void AddVariable(const hephaestus::InputParameters);
  void Init(mfem::ParMesh &pmesh);

  std::vector<hephaestus::InputParameters> variable_params;
  mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> fespaces;
  mfem::NamedFieldsMap<mfem::ParGridFunction> gfs;
};

} // namespace hephaestus
