#include "equation_system_operator.hpp"

namespace hephaestus {

void EquationSystemOperator::SetGridFunctions() {
  local_test_vars = populateVectorFromNamedFieldsMap<mfem::ParGridFunction>(
      _gridfunctions, state_var_names);

  // Set operator size and block structure
  block_trueOffsets.SetSize(local_test_vars.size() + 1);
  block_trueOffsets[0] = 0;
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    block_trueOffsets[ind + 1] =
        local_test_vars.at(ind)->ParFESpace()->TrueVSize();
  }
  block_trueOffsets.PartialSum();

  true_offsets.SetSize(local_test_vars.size() + 1);
  true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    true_offsets[ind + 1] = local_test_vars.at(ind)->ParFESpace()->GetVSize();
  }
  true_offsets.PartialSum();

  this->height = true_offsets[local_test_vars.size()];
  this->width = true_offsets[local_test_vars.size()];
  trueX.Update(block_trueOffsets);
  trueRhs.Update(block_trueOffsets);

  // Populate vector of active auxiliary gridfunctions
  active_aux_var_names.resize(0);
  for (auto &aux_var_name : aux_var_names) {
    if (_gridfunctions.Has(aux_var_name)) {
      active_aux_var_names.push_back(aux_var_name);
    }
  }
};

void EquationSystemOperator::Init(mfem::Vector &X) {
  // Define material property coefficients
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    local_test_vars.at(ind)->MakeRef(local_test_vars.at(ind)->ParFESpace(),
                                     const_cast<mfem::Vector &>(X),
                                     true_offsets[ind]);
  }
}

void EquationSystemOperator::Solve(mfem::Vector &X) {}

} // namespace hephaestus
