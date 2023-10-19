#include "time_domain_equation_system_operator.hpp"

namespace hephaestus {

std::string GetTimeDerivativeName(const std::string &name) {
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

void TimeDomainEquationSystemOperator::SetGridFunctions() {
  state_var_names = _equation_system->var_names;
  local_test_vars = populateVectorFromNamedFieldsMap<mfem::ParGridFunction>(
      _gridfunctions, _equation_system->var_names);
  local_trial_vars = populateVectorFromNamedFieldsMap<mfem::ParGridFunction>(
      _gridfunctions, _equation_system->var_time_derivative_names);

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

void TimeDomainEquationSystemOperator::Init(mfem::Vector &X) {
  // Define material property coefficients
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    local_test_vars.at(ind)->MakeRef(local_test_vars.at(ind)->ParFESpace(),
                                     const_cast<mfem::Vector &>(X),
                                     true_offsets[ind]);
    *(local_test_vars.at(ind)) = 0.0;
    *(local_trial_vars.at(ind)) = 0.0;
  }

  _equation_system->buildEquationSystem(_bc_map, _sources);
};

void TimeDomainEquationSystemOperator::ImplicitSolve(const double dt,
                                                     const mfem::Vector &X,
                                                     mfem::Vector &dX_dt) {
  dX_dt = 0.0;
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    local_test_vars.at(ind)->MakeRef(local_test_vars.at(ind)->ParFESpace(),
                                     const_cast<mfem::Vector &>(X),
                                     true_offsets[ind]);
    local_trial_vars.at(ind)->MakeRef(local_trial_vars.at(ind)->ParFESpace(),
                                      dX_dt, true_offsets[ind]);
  }
  _coefficients.SetTime(this->GetTime());
  _equation_system->setTimeStep(dt);
  _equation_system->updateEquationSystem(_bc_map, _sources);

  _equation_system->FormLinearSystem(blockA, trueX, trueRhs);
  if (solver != NULL) {
    delete solver;
  }
  solver = new hephaestus::DefaultGMRESSolver(
      _solver_options, *blockA.As<mfem::HypreParMatrix>());
  solver->Mult(trueRhs, trueX);
  _equation_system->RecoverFEMSolution(trueX, _gridfunctions);
}

void TimeDomainEquationSystemOperator::SetEquationSystem(
    hephaestus::TimeDependentEquationSystem *equation_system) {
  _equation_system = equation_system;
}

} // namespace hephaestus
