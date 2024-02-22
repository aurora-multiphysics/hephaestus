#include "steady_state_problem_builder.hpp"

namespace hephaestus
{

void
SteadyStateProblemBuilder::ConstructEquationSystem()
{
  hephaestus::InputParameters params;
  auto equation_system = std::make_unique<hephaestus::EquationSystem>(params);

  _problem->GetOperator()->SetEquationSystem(std::move(equation_system));
}

void
SteadyStateProblemBuilder::SetOperatorGridFunctions()
{
  _problem->GetOperator()->SetGridFunctions();
}

void
SteadyStateProblemBuilder::InitializeKernels()
{
  if (_problem->HasEquationSystem())
  {
    _problem->GetEquationSystem()->Init(
        _problem->_gridfunctions, _problem->_fespaces, _problem->_bc_map, _problem->_coefficients);
  }

  _problem->_preprocessors.Init(_problem->_gridfunctions, _problem->_coefficients);
  _problem->_sources.Init(
      _problem->_gridfunctions, _problem->_fespaces, _problem->_bc_map, _problem->_coefficients);
}

void
SteadyStateProblemBuilder::ConstructOperator()
{
  _problem->SetOperator(std::make_unique<hephaestus::ProblemOperator>(*_problem));
}

void
SteadyStateProblemBuilder::ConstructState()
{
  _problem->_f =
      std::make_unique<mfem::BlockVector>(_problem->GetOperator()->_true_offsets); // Vector of dofs
  _problem->GetOperator()->Init(*(_problem->_f)); // Set up initial conditions
}
} // namespace hephaestus
