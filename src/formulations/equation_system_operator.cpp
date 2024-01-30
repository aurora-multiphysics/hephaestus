#include "equation_system_operator.hpp"

namespace hephaestus
{

void
EquationSystemOperator::SetGridFunctions()
{
  _local_test_vars = _gridfunctions.Get(_state_var_names);

  // Set operator size and block structure
  _block_true_offsets.SetSize(_local_test_vars.size() + 1);
  _block_true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < _local_test_vars.size(); ++ind)
  {
    _block_true_offsets[ind + 1] = _local_test_vars.at(ind)->ParFESpace()->TrueVSize();
  }
  _block_true_offsets.PartialSum();

  _true_offsets.SetSize(_local_test_vars.size() + 1);
  _true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < _local_test_vars.size(); ++ind)
  {
    _true_offsets[ind + 1] = _local_test_vars.at(ind)->ParFESpace()->GetVSize();
  }
  _true_offsets.PartialSum();

  height = _true_offsets[_local_test_vars.size()];
  width = _true_offsets[_local_test_vars.size()];
  _true_x.Update(_block_true_offsets);
  _true_rhs.Update(_block_true_offsets);

  // Populate vector of active auxiliary gridfunctions
  _active_aux_var_names.resize(0);
  for (auto & aux_var_name : _aux_var_names)
  {
    if (_gridfunctions.Has(aux_var_name))
    {
      _active_aux_var_names.push_back(aux_var_name);
    }
  }
};

void
EquationSystemOperator::SetSolver(const mfem::HypreParMatrix & M)
{
  auto solve_type = _solver_options.GetParam<std::string>("Solver");

  if (solve_type == "PCG")
    _solver = std::make_unique<DefaultH1PCGSolver>(_solver_options, M);
  else if (solve_type == "GMRES")
    _solver = std::make_unique<DefaultGMRESSolver>(_solver_options, M);
  else if (solve_type == "Jacobi")
    _solver = std::make_unique<DefaultJacobiPCGSolver>(_solver_options, M);
  else
    mfem::mfem_error("H1 Solver type not recognised! Please set Solver field in solver_options "
                     "to one of the following values: PCG, GMRES, Jacobi");
}

void
EquationSystemOperator::SetSolver(const mfem::HypreParMatrix & M,
                                  mfem::ParFiniteElementSpace * edge_fespace)
{
  auto solve_type = _solver_options.GetParam<std::string>("Solver");

  if (solve_type == "HCurl_PCG")
    _solver = std::make_unique<DefaultHCurlPCGSolver>(_solver_options, M, edge_fespace);
  else if (solve_type == "HCurl_FGMRES")
    _solver = std::make_unique<DefaultHCurlFGMRESSolver>(_solver_options, M, edge_fespace);
  else
    mfem::mfem_error("HCurl Solver type not recognised! Please set Solver field in "
                     "solver_options to one of the following values: HCurl_PCG, HCurl_FGMRES");
}

void
EquationSystemOperator::Init(mfem::Vector & X)
{
  // Define material property coefficients
  for (unsigned int ind = 0; ind < _local_test_vars.size(); ++ind)
  {
    _local_test_vars.at(ind)->MakeRef(
        _local_test_vars.at(ind)->ParFESpace(), const_cast<mfem::Vector &>(X), _true_offsets[ind]);
  }
}

void
EquationSystemOperator::Solve(mfem::Vector & X)
{
}

} // namespace hephaestus
