#include "variables.hpp"

namespace hephaestus {

void FESpaces::StoreInput(const hephaestus::InputParameters var_params) {
  input_params.push_back(var_params);
}

void FESpaces::Init(mfem::ParMesh &pmesh) {
  for (const auto &params : input_params) {
    mfem::ParFiniteElementSpace *parfespace =
        hephaestus::Factory::createParFESpace(params, pmesh);
    Register(params.GetParam<std::string>("FESpaceName"), parfespace, true);
  }
}

void GridFunctions::Init(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces) {

  for (const auto &params : input_params) {
    mfem::ParFiniteElementSpace *fespace(
        fespaces.Get(params.GetParam<std::string>("FESpaceName")));

    mfem::ParGridFunction *gridfunc = new mfem::ParGridFunction(
        fespaces.Get(params.GetParam<std::string>("FESpaceName")));

    *gridfunc = 0.0;

    Register(params.GetParam<std::string>("VariableName"), gridfunc, true);
  }
}

void GridFunctions::StoreInput(const hephaestus::InputParameters var_params) {
  input_params.push_back(var_params);
}

} // namespace hephaestus
