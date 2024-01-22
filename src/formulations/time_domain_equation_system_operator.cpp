#include "time_domain_equation_system_operator.hpp"

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
TimeDomainEquationSystemOperator::SetGridFunctions()
{
  _state_var_names = _equation_system->_var_names;
  _local_test_vars = _gridfunctions.Get(_equation_system->_var_names);
  _local_trial_vars = _gridfunctions.Get(_equation_system->_var_time_derivative_names);

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
TimeDomainEquationSystemOperator::Init(mfem::Vector & X)
{
  // Define material property coefficients
  for (unsigned int ind = 0; ind < _local_test_vars.size(); ++ind)
  {
    _local_test_vars.at(ind)->MakeRef(
        _local_test_vars.at(ind)->ParFESpace(), const_cast<mfem::Vector &>(X), _true_offsets[ind]);
    *(_local_test_vars.at(ind)) = 0.0;
    *(_local_trial_vars.at(ind)) = 0.0;
  }

  _equation_system->BuildEquationSystem(_bc_map, _sources);
};

void
TimeDomainEquationSystemOperator::ImplicitSolve(const double dt,
                                                const mfem::Vector & X,
                                                mfem::Vector & dX_dt)
{
  dX_dt = 0.0;
  for (unsigned int ind = 0; ind < _local_test_vars.size(); ++ind)
  {
    _local_test_vars.at(ind)->MakeRef(
        _local_test_vars.at(ind)->ParFESpace(), const_cast<mfem::Vector &>(X), _true_offsets[ind]);
    _local_trial_vars.at(ind)->MakeRef(
        _local_trial_vars.at(ind)->ParFESpace(), dX_dt, _true_offsets[ind]);
  }
  _coefficients.SetTime(GetTime());
  _equation_system->SetTimeStep(dt);
  _equation_system->UpdateEquationSystem(_bc_map, _sources);

  _equation_system->FormLinearSystem(_block_a, _true_x, _true_rhs);

  _solver = std::make_unique<hephaestus::DefaultGMRESSolver>(_solver_options,
                                                             *_block_a.As<mfem::HypreParMatrix>());

  _solver->Mult(_true_rhs, _true_x);
  _equation_system->RecoverFEMSolution(_true_x, _gridfunctions);
}

void
TimeDomainEquationSystemOperator::SetEquationSystem(
    hephaestus::TimeDependentEquationSystem * equation_system)
{
  _equation_system = equation_system;
}

} // namespace hephaestus
