#include "problem_operator.hpp"

namespace hephaestus {

void ProblemOperator::SetGridFunctions() {
  trial_variables = populateVectorFromNamedFieldsMap<mfem::ParGridFunction>(
      _gridfunctions, trial_var_names);

  // Set operator size and block structure
  block_trueOffsets.SetSize(trial_variables.size() + 1);
  block_trueOffsets[0] = 0;
  for (unsigned int ind = 0; ind < trial_variables.size(); ++ind) {
    block_trueOffsets[ind + 1] =
        trial_variables.at(ind)->ParFESpace()->TrueVSize();
  }
  block_trueOffsets.PartialSum();

  true_offsets.SetSize(trial_variables.size() + 1);
  true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < trial_variables.size(); ++ind) {
    true_offsets[ind + 1] = trial_variables.at(ind)->ParFESpace()->GetVSize();
  }
  true_offsets.PartialSum();

  height = true_offsets[trial_variables.size()];
  width = true_offsets[trial_variables.size()];
  trueX.Update(block_trueOffsets);
  trueRhs.Update(block_trueOffsets);
};

void ProblemOperator::Init(mfem::Vector &X) {
  // Define material property coefficients
  for (unsigned int ind = 0; ind < trial_variables.size(); ++ind) {
    trial_variables.at(ind)->MakeRef(trial_variables.at(ind)->ParFESpace(),
                                     const_cast<mfem::Vector &>(X),
                                     true_offsets[ind]);
  }
}

} // namespace hephaestus
