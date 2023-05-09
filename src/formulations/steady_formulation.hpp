#pragma once
#include "../common/pfem_extras.hpp"
#include "auxkernels.hpp"
#include "equation_system.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {
// Specifies output interfaces of a time-domain EM formulation.
//   Curl mu^{-1} Curl E - omega^2 epsilon E + i omega sigma E = - i omega J
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

//
// Specifies output interfaces of a time-domain EM formulation.
class SteadyFormulation {
  // std::vector<mfem::ParGridFunction *> local_trial_vars, local_test_vars;
  std::string frequency_coef_name, permittivity_coef_name,
      reluctivity_coef_name, conductivity_coef_name, h_curl_var_name;

public:
  hephaestus::EquationSystem *equation_system;
  hephaestus::FrequencyDomainOperator *fd_operator;
  mfem::ConstantCoefficient oneCoef;
  mfem::ConstantCoefficient *freqCoef;

  SteadyFormulation();

  virtual hephaestus::FrequencyDomainOperator *CreateFrequencyDomainOperator(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options) {
    fd_operator = new hephaestus::FrequencyDomainOperator(
        pmesh, fespaces, variables, bc_map, domain_properties, sources,
        solver_options);
    return fd_operator;
  };

  virtual void RegisterMissingVariables(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables){};

  virtual void
  RegisterAuxKernels(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                     hephaestus::AuxKernels &auxkernels){};

  virtual void
  RegisterCoefficients(hephaestus::DomainProperties &domain_properties){};
};
} // namespace hephaestus
