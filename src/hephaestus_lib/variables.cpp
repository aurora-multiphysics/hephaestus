#include "variables.hpp"

namespace hephaestus {

void Variables::AddVariable(const hephaestus::InputParameters var_params) {
  variable_params.push_back(var_params);
}

void Variables::Init(mfem::ParMesh &pmesh) {
  mfem::Vector zero_vec(3);
  zero_vec = 0.0;
  mfem::VectorConstantCoefficient Zero_vec(zero_vec);
  mfem::ConstantCoefficient Zero(0.0);

  for (const auto &params : variable_params) {
    mfem::ParFiniteElementSpace *parfespace =
        hephaestus::Factory::createParFESpace(params, pmesh);

    fespaces.Register(params.GetParam<std::string>("FESpaceName"), parfespace,
                      true);

    mfem::ParGridFunction *gridfunc = new mfem::ParGridFunction(
        fespaces.Get(params.GetParam<std::string>("FESpaceName")));

    std::string fespacetype = params.GetParam<std::string>("FESpaceType");
    if (fespacetype == std::string("H1") || fespacetype == std::string("L2")) {
      gridfunc->ProjectCoefficient(Zero);
    } else {
      gridfunc->ProjectCoefficient(Zero_vec);
    }

    gfs.Register(params.GetParam<std::string>("VariableName"), gridfunc, true);
  }
}

} // namespace hephaestus
