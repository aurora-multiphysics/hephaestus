#pragma once
#include "../common/pfem_extras.hpp"
#include "factory.hpp"
#include "inputs.hpp"
#include "mfem.hpp"

namespace hephaestus {

class FESpaces : public mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> {
public:
  void StoreInput(const hephaestus::InputParameters);
  void Init(mfem::ParMesh &pmesh);

  std::vector<hephaestus::InputParameters> input_params;
};

class GridFunctions : public mfem::NamedFieldsMap<mfem::ParGridFunction> {
public:
  void StoreInput(const hephaestus::InputParameters);
  void Init(mfem::ParMesh &pmesh,
            mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces);

  std::vector<hephaestus::InputParameters> input_params;
};

} // namespace hephaestus
