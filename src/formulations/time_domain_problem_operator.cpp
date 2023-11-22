#include "time_domain_problem_operator.hpp"

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

void TimeDomainProblemOperator::SetGridFunctions() {
  trial_var_names = _equation_system->trial_var_names;
  trial_variables = populateVectorFromNamedFieldsMap<mfem::ParGridFunction>(
      _gridfunctions, _equation_system->trial_var_names);
  trial_variable_time_derivatives =
      populateVectorFromNamedFieldsMap<mfem::ParGridFunction>(
          _gridfunctions, _equation_system->trial_var_time_derivative_names);

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

  this->height = true_offsets[trial_variables.size()];
  this->width = true_offsets[trial_variables.size()];
  trueX.Update(block_trueOffsets);
  trueRhs.Update(block_trueOffsets);
};

void TimeDomainProblemOperator::Init(mfem::Vector &X) {
  // Define material property coefficients
  for (unsigned int ind = 0; ind < trial_variables.size(); ++ind) {
    trial_variables.at(ind)->MakeRef(trial_variables.at(ind)->ParFESpace(),
                                     const_cast<mfem::Vector &>(X),
                                     true_offsets[ind]);
    *(trial_variables.at(ind)) = 0.0;
    *(trial_variable_time_derivatives.at(ind)) = 0.0;
  }

  _equation_system->buildEquationSystem(_bc_map, _sources);
};

void TimeDomainProblemOperator::ImplicitSolve(const double dt,
                                              const mfem::Vector &X,
                                              mfem::Vector &dX_dt) {
  dX_dt = 0.0;
  for (unsigned int ind = 0; ind < trial_variables.size(); ++ind) {
    trial_variables.at(ind)->MakeRef(trial_variables.at(ind)->ParFESpace(),
                                     const_cast<mfem::Vector &>(X),
                                     true_offsets[ind]);
    trial_variable_time_derivatives.at(ind)->MakeRef(
        trial_variable_time_derivatives.at(ind)->ParFESpace(), dX_dt,
        true_offsets[ind]);
  }
  _coefficients.SetTime(this->GetTime());
  buildEquationSystemOperator(dt);
  buildJacobianSolver();
  getJacobianSolver()->Mult(trueRhs, trueX);
  _equation_system->RecoverFEMSolution(trueX, _gridfunctions);
}
void TimeDomainProblemOperator::buildEquationSystemOperator(double dt) {
  _equation_system->setTimeStep(dt);
  _equation_system->updateEquationSystem(_bc_map, _sources);
  _equation_system->FormLinearSystem(_equation_system_operator, trueX, trueRhs);
}
void TimeDomainProblemOperator::buildJacobianSolver() {
  setJacobianSolver(new hephaestus::DefaultGMRESSolver(
      _solver_options, *_equation_system_operator.As<mfem::HypreParMatrix>()));
}

void TimeDomainProblemOperator::SetEquationSystem(
    hephaestus::TimeDependentEquationSystem *equation_system) {
  _equation_system = equation_system;
}

} // namespace hephaestus
