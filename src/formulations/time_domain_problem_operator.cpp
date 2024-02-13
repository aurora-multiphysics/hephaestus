#include "time_domain_problem_operator.hpp"

namespace hephaestus
{

std::string
GetTimeDerivativeName(const std::string & name)
{
  return std::string("d") + name + std::string("_dt");
}

std::vector<std::string>
GetTimeDerivativeNames(std::vector<std::string> gridfunction_names)
{
  std::vector<std::string> time_derivative_names;
  for (auto & gridfunction_name : gridfunction_names)
  {
    time_derivative_names.push_back(GetTimeDerivativeName(gridfunction_name));
  }
  return time_derivative_names;
}

void
TimeDomainProblemOperator::SetGridFunctions()
{
  _trial_var_names = _equation_system->_trial_var_names;
  _trial_variables = _problem._gridfunctions.Get(_equation_system->_trial_var_names);
  _trial_variable_time_derivatives =
      _problem._gridfunctions.Get(_equation_system->_trial_var_time_derivative_names);

  // Set operator size and block structure
  _block_true_offsets.SetSize(_trial_variables.size() + 1);
  _block_true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < _trial_variables.size(); ++ind)
  {
    _block_true_offsets[ind + 1] = _trial_variables.at(ind)->ParFESpace()->TrueVSize();
  }
  _block_true_offsets.PartialSum();

  _true_offsets.SetSize(_trial_variables.size() + 1);
  _true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < _trial_variables.size(); ++ind)
  {
    _true_offsets[ind + 1] = _trial_variables.at(ind)->ParFESpace()->GetVSize();
  }
  _true_offsets.PartialSum();

  height = _true_offsets[_trial_variables.size()];
  width = _true_offsets[_trial_variables.size()];
  _true_x.Update(_block_true_offsets);
  _true_rhs.Update(_block_true_offsets);
};

void
TimeDomainProblemOperator::Init(mfem::Vector & X)
{
  // Define material property coefficients
  for (unsigned int ind = 0; ind < _trial_variables.size(); ++ind)
  {
    _trial_variables.at(ind)->MakeRef(
        _trial_variables.at(ind)->ParFESpace(), const_cast<mfem::Vector &>(X), _true_offsets[ind]);
    *(_trial_variables.at(ind)) = 0.0;
    *(_trial_variable_time_derivatives.at(ind)) = 0.0;
  }

  _equation_system->BuildEquationSystem(_problem._bc_map, _problem._sources);
};

void
TimeDomainProblemOperator::ImplicitSolve(const double dt,
                                         const mfem::Vector & X,
                                         mfem::Vector & dX_dt)
{
  dX_dt = 0.0;
  for (unsigned int ind = 0; ind < _trial_variables.size(); ++ind)
  {
    _trial_variables.at(ind)->MakeRef(
        _trial_variables.at(ind)->ParFESpace(), const_cast<mfem::Vector &>(X), _true_offsets[ind]);
    _trial_variable_time_derivatives.at(ind)->MakeRef(
        _trial_variable_time_derivatives.at(ind)->ParFESpace(), dX_dt, _true_offsets[ind]);
  }
  _problem._coefficients.SetTime(GetTime());
  BuildEquationSystemOperator(dt);

  _problem._nonlinear_solver->SetSolver(*_problem._jacobian_solver);
  _problem._nonlinear_solver->SetOperator(*_equation_system);
  _problem._nonlinear_solver->Mult(_true_rhs, _true_x);

  _equation_system->RecoverFEMSolution(_true_x, _problem._gridfunctions);
}
void
TimeDomainProblemOperator::BuildEquationSystemOperator(double dt)
{
  _equation_system->SetTimeStep(dt);
  _equation_system->UpdateEquationSystem(_problem._bc_map, _problem._sources);
  _equation_system->BuildJacobian(_true_x, _true_rhs);
}

} // namespace hephaestus