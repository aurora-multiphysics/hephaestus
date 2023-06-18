#pragma once
#include "formulation_base.hpp"

namespace hephaestus {
class FrequencyDomainOperator : public mfem::Operator {
public:
  FrequencyDomainOperator(
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

  ~FrequencyDomainOperator(){};

  virtual void SetVariables();
  virtual void Init(mfem::Vector &X);
  virtual void Solve(mfem::Vector &X);
  void Mult(const mfem::Vector &x, mfem::Vector &y) const override{};

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

  int myid_;
  int num_procs_;
  mfem::ParMesh *pmesh_;
  mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &_fespaces;
  mfem::NamedFieldsMap<mfem::ParGridFunction> &_variables;
  hephaestus::BCMap &_bc_map;
  hephaestus::Sources &_sources;
  hephaestus::DomainProperties &_domain_properties;
  hephaestus::InputParameters &_solver_options;

  mfem::OperatorHandle blockA;
  mfem::BlockVector trueX, trueRhs;
};

} // namespace hephaestus
