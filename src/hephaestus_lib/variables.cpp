#include "variables.hpp"

namespace hephaestus {

void Variables::AddVariable(const hephaestus::InputParameters var_params) {
  variable_params.push_back(var_params);
}

void Variables::Init(mfem::ParMesh &pmesh) {
  for (const auto &params : variable_params) {
    mfem::ParFiniteElementSpace *parfespace =
        hephaestus::Factory::createParFESpace(params, pmesh);
    fespaces.Register(params.GetParam<std::string>("FESpaceName"), parfespace,
                      true);
    gfs.Register(params.GetParam<std::string>("VariableName"),
                 new mfem::ParGridFunction(parfespace), true);
  }
}

} // namespace hephaestus
