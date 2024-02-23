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
TimeDomainEquationSystemProblemOperator::SetGridFunctions()
{
  _trial_var_names = GetEquationSystem()->_trial_var_names;
  _trial_variables = _problem._gridfunctions.Get(GetEquationSystem()->_trial_var_names);
  _trial_variable_time_derivatives =
      _problem._gridfunctions.Get(GetEquationSystem()->_trial_var_time_derivative_names);

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
}

void
TimeDomainEquationSystemProblemOperator::Init(mfem::Vector & X)
{
  // Define material property coefficients
  for (unsigned int ind = 0; ind < _trial_variables.size(); ++ind)
  {
    _trial_variables.at(ind)->MakeRef(
        _trial_variables.at(ind)->ParFESpace(), const_cast<mfem::Vector &>(X), _true_offsets[ind]);
    *(_trial_variables.at(ind)) = 0.0;
    *(_trial_variable_time_derivatives.at(ind)) = 0.0;
  }

  GetEquationSystem()->BuildEquationSystem(_problem._bc_map, _problem._sources);
}

void
TimeDomainEquationSystemProblemOperator::ImplicitSolve(const double dt,
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
  _problem._nonlinear_solver->SetOperator(*GetEquationSystem());
  _problem._nonlinear_solver->Mult(_true_rhs, _true_x);

  GetEquationSystem()->RecoverFEMSolution(_true_x, _problem._gridfunctions);
}

void
TimeDomainEquationSystemProblemOperator::BuildEquationSystemOperator(double dt)
{
  GetEquationSystem()->SetTimeStep(dt);
  GetEquationSystem()->UpdateEquationSystem(_problem._bc_map, _problem._sources);
  GetEquationSystem()->BuildJacobian(_true_x, _true_rhs);
}

} // namespace hephaestus