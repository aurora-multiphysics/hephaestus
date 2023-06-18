#pragma once
#include "formulation_base.hpp"

namespace hephaestus {

std::string GetTimeDerivativeName(std::string name);

std::vector<std::string>
GetTimeDerivativeNames(std::vector<std::string> gridfunction_names);

// Specifies output interfaces of a time-domain EM formulation.
class TimeDomainEquationSystemOperator : public mfem::TimeDependentOperator {
public:
  TimeDomainEquationSystemOperator(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
      : myid_(0), num_procs_(1), pmesh_(&pmesh), _fespaces(fespaces),
        _variables(variables), _bc_map(bc_map), _sources(sources),
        _domain_properties(domain_properties),
        _solver_options(solver_options){};

  ~TimeDomainEquationSystemOperator(){};

  virtual void SetVariables();
  virtual void Init(mfem::Vector &X);
  virtual void ImplicitSolve(const double dt, const mfem::Vector &X,
                             mfem::Vector &dX_dt) override;
  void
  SetEquationSystem(hephaestus::TimeDependentEquationSystem *equation_system);

  mfem::Array<int> true_offsets, block_trueOffsets;
  // Vector of names of state variables used in formulation, ordered by
  // appearance in block vector during solve.
  std::vector<std::string> state_var_names;
  // Vector of names of recognised auxiliary variables that can be calculated
  // from formulation,
  std::vector<std::string> aux_var_names;
  // Vector of names of active auxiliary variables that are being calculated
  // in formulation,
  std::vector<std::string> active_aux_var_names;

  std::vector<mfem::ParGridFunction *> local_trial_vars, local_test_vars;

  hephaestus::TimeDependentEquationSystem *_equation_system;

  int myid_;
  int num_procs_;
  mfem::ParMesh *pmesh_;
  mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &_fespaces;
  mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables;
  hephaestus::BCMap &_bc_map;
  hephaestus::Sources &_sources;
  hephaestus::DomainProperties &_domain_properties;
  hephaestus::InputParameters &_solver_options;

  mutable hephaestus::DefaultGMRESSolver *solver = NULL;
  mutable hephaestus::DefaultHCurlPCGSolver *a1_solver = NULL;

  mfem::OperatorHandle blockA;
  mfem::BlockVector trueX, trueRhs;
};

} // namespace hephaestus
