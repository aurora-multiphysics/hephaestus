#include "time_domain_problem_builder.hpp"

namespace hephaestus
{

std::vector<mfem::ParGridFunction *>
TimeDomainProblemBuilder::RegisterTimeDerivatives(std::vector<std::string> gridfunction_names,
                                                  hephaestus::GridFunctions & gridfunctions)
{
  std::vector<mfem::ParGridFunction *> time_derivatives;

  for (auto & gridfunction_name : gridfunction_names)
  {
    gridfunctions.Register(GetTimeDerivativeName(gridfunction_name),
                           std::make_shared<mfem::ParGridFunction>(
                               gridfunctions.Get(gridfunction_name)->ParFESpace()));

    time_derivatives.push_back(gridfunctions.Get(GetTimeDerivativeName(gridfunction_name)));
  }

  return time_derivatives;
}

void
TimeDomainProblemBuilder::RegisterGridFunctions()
{
  std::vector<std::string> gridfunction_names;
  for (auto const & [name, gf] : GetProblem()->_gridfunctions)
  {
    gridfunction_names.push_back(name);
  }
  RegisterTimeDerivatives(gridfunction_names, GetProblem()->_gridfunctions);
}

void
TimeDomainProblemBuilder::SetOperatorGridFunctions()
{
  GetProblem()->GetOperator()->SetGridFunctions();
}

void
TimeDomainProblemBuilder::ConstructOperator()
{
  GetProblem()->ConstructOperator();
}

void
TimeDomainProblemBuilder::ConstructState()
{
  ProblemBuilder::ConstructState();
  GetProblem()->GetOperator()->SetTime(0.0);
}

void
TimeDomainProblemBuilder::ConstructTimestepper()
{
  GetProblem()->_ode_solver = std::make_unique<mfem::BackwardEulerSolver>();
  GetProblem()->_ode_solver->Init(*(GetProblem()->GetOperator()));
}

} // namespace hephaestus
