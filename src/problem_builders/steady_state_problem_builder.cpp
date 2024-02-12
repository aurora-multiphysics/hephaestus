#include "steady_state_problem_builder.hpp"

namespace hephaestus
{

void
SteadyStateProblemBuilder::InitializeKernels()
{
  _problem->_preprocessors.Init(_problem->_gridfunctions, _problem->_coefficients);
  _problem->_sources.Init(
      _problem->_gridfunctions, _problem->_fespaces, _problem->_bc_map, _problem->_coefficients);
}

void
SteadyStateProblemBuilder::ConstructOperator()
{
  _problem->_ss_operator = std::make_unique<hephaestus::ProblemOperator>(*_problem);
  _problem->_ss_operator->SetGridFunctions();
}

void
SteadyStateProblemBuilder::ConstructState()
{
  _problem->_f =
      std::make_unique<mfem::BlockVector>(_problem->_ss_operator->_true_offsets); // Vector of dofs
  _problem->_ss_operator->Init(*(_problem->_f)); // Set up initial conditions
}
} // namespace hephaestus
