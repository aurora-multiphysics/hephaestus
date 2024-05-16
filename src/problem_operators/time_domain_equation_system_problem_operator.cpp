#include "time_domain_equation_system_problem_operator.hpp"

namespace hephaestus
{

void
TimeDomainEquationSystemProblemOperator::SetTrialVariableNames()
{
  _trial_var_names = GetEquationSystem()->_trial_var_names;
}

void
TimeDomainEquationSystemProblemOperator::SetTrialVariables()
{
  TimeDomainProblemOperator::SetTrialVariables();

  _trial_variable_time_derivatives =
      _problem._gridfunctions.Get(GetEquationSystem()->_trial_var_time_derivative_names);

  // Define material property coefficients
  for (size_t i = 0; i < _trial_variables.size(); ++i)
  {
    *(_trial_variable_time_derivatives.at(i)) = 0.0;
  }
}

void
TimeDomainEquationSystemProblemOperator::Init()
{
  GetEquationSystem()->Init(_problem._gridfunctions,
                            _problem._fespaces,
                            _problem._bc_map,
                            _problem._coefficients,
                            _problem._sources);

  TimeDomainProblemOperator::Init();
}

void
TimeDomainEquationSystemProblemOperator::Update()
{
  GetEquationSystem()->Update(_problem._bc_map, _problem._sources);

  TimeDomainProblemOperator::Update();

  //
  // TODO: - we need to update the size of the jacobian_solver here after the parent class' Update
  // method is called which ensures that we've updated the _true_x, _true_rhs. This is a hacky first
  // attempt to avoid a segfault with the test_mesh_updates. This is very bad code and basically
  // copies methods from the problem builder used. We need to generalize this.

  GetEquationSystem()->BuildJacobian(_true_x, _true_rhs);

  // Rebuild the jacobian preconditioner.
  auto first_pfespace = GetEquationSystem()->_test_pfespaces.at(0);

  auto precond = std::make_shared<mfem::HypreAMS>(first_pfespace);

  precond->SetSingularProblem();
  precond->SetPrintLevel(-1);

  _problem._jacobian_preconditioner = precond;

  // Set new preconditioner.
  std::dynamic_pointer_cast<mfem::HyprePCG>(_problem._jacobian_solver)
      ->SetPreconditioner(
          *std::dynamic_pointer_cast<mfem::HypreSolver>(_problem._jacobian_preconditioner));

  // Set Jacobian matrix.
  auto * matrix = GetEquationSystem()->JacobianOperatorHandle().As<mfem::HypreParMatrix>();
  _problem._jacobian_solver->SetOperator(*matrix);
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

  // TODO: - We have a problem here on the second ImplicitSolve after an update. This is due to the
  // jacobian solver being of the incorrect size. Fails mysteriously somewhere in the Mult method.

  _problem._nonlinear_solver->SetSolver(*_problem._jacobian_solver);
  _problem._nonlinear_solver->SetOperator(*GetEquationSystem());
  _problem._nonlinear_solver->Mult(_true_rhs, _true_x);

  GetEquationSystem()->RecoverFEMSolution(_true_x, _problem._gridfunctions);
}

void
TimeDomainEquationSystemProblemOperator::BuildEquationSystemOperator(double dt)
{
  GetEquationSystem()->SetTimeStep(dt);
  GetEquationSystem()->BuildEquationSystem(_problem._bc_map, _problem._sources);
  GetEquationSystem()->BuildJacobian(_true_x, _true_rhs);
}

} // namespace hephaestus