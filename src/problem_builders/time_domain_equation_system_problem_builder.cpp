#include "time_domain_equation_system_problem_builder.hpp"

namespace hephaestus
{

void
TimeDomainEquationSystemProblemBuilder::InitializeEquationSystem()
{
  GetEquationSystem()->Init(GetProblem()->_gridfunctions,
                            GetProblem()->_fespaces,
                            GetProblem()->_bc_map,
                            GetProblem()->_coefficients);
}

} // namespace hephaestus
