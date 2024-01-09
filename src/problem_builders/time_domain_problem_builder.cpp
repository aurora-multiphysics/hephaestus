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
    gridfunctions.Register(
        GetTimeDerivativeName(gridfunction_name),
        new mfem::ParGridFunction(gridfunctions.Get(gridfunction_name)->ParFESpace()),
        true);
    time_derivatives.push_back(gridfunctions.Get(GetTimeDerivativeName(gridfunction_name)));
  }
  return time_derivatives;
}

void
TimeDomainProblemBuilder::RegisterGridFunctions()
{
  std::vector<std::string> gridfunction_names;
  for (auto const & [name, gf] : problem->gridfunctions)
  {
    gridfunction_names.push_back(name);
  }
  RegisterTimeDerivatives(gridfunction_names, problem->gridfunctions);
}

void
TimeDomainProblemBuilder::ConstructEquationSystem()
{
  hephaestus::InputParameters params;
  problem->td_equation_system = std::make_unique<hephaestus::TimeDependentEquationSystem>(params);
}

void
TimeDomainProblemBuilder::InitializeKernels()
{
  problem->td_equation_system->Init(
      problem->gridfunctions, problem->fespaces, problem->bc_map, problem->coefficients);

  problem->preprocessors.Init(problem->gridfunctions, problem->coefficients);
  problem->sources.Init(
      problem->gridfunctions, problem->fespaces, problem->bc_map, problem->coefficients);
}

void
TimeDomainProblemBuilder::ConstructOperator()
{
  problem->td_operator =
      std::make_unique<hephaestus::TimeDomainEquationSystemOperator>(*(problem->pmesh),
                                                                     problem->fespaces,
                                                                     problem->gridfunctions,
                                                                     problem->bc_map,
                                                                     problem->coefficients,
                                                                     problem->sources,
                                                                     problem->solver_options);
  problem->td_operator->SetEquationSystem(problem->td_equation_system.get());
  problem->td_operator->SetGridFunctions();
}

void
TimeDomainProblemBuilder::ConstructState()
{
  problem->F = new mfem::BlockVector(problem->td_operator->true_offsets); // Vector of dofs
  *(problem->F) = 0.0;                                                    // give initial value
  problem->td_operator->Init(*(problem->F)); // Set up initial conditions
  problem->td_operator->SetTime(0.0);
}

void
TimeDomainProblemBuilder::ConstructSolver()
{
  problem->ode_solver = new mfem::BackwardEulerSolver;
  problem->ode_solver->Init(*(problem->td_operator));
}

} // namespace hephaestus
