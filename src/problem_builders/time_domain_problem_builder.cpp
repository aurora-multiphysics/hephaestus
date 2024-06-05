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
TimeDomainProblemBuilder::ConstructOperator()
{
  InputParameters params;
  params.Set("Problem", GetBaseProblem());

  auto problem_operator = std::make_unique<hephaestus::TimeDomainProblemOperator>(params);
  GetProblem()->SetOperator(std::move(problem_operator));
}

void
TimeDomainProblemBuilder::InitializeOperator()
{
  ProblemBuilder::InitializeOperator();
  GetProblem()->GetOperator()->SetTime(0.0);
}

} // namespace hephaestus
