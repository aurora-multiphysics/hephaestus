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

TimeDomainProblemOperator::TimeDomainProblemOperator(const InputParameters & params)
  : ProblemOperatorBase(params)
{
}

void
TimeDomainProblemOperator::ConstructTimestepper()
{
  _ode_solver = std::make_unique<mfem::BackwardEulerSolver>();
  _ode_solver->Init(*this);
}

void
TimeDomainProblemOperator::Init()
{
  ProblemOperatorBase::Init();
  ConstructTimestepper();
}

void
TimeDomainProblemOperator::Update()
{
  ProblemOperatorBase::Update();

  // The dimensions of the problem operator have now changed. We must call the
  // ODE solver's Init method to resize its internal vector.
  _ode_solver->Init(*this);
}

} // namespace hephaestus