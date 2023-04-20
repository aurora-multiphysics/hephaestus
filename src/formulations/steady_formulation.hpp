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
class HertzOperator : public mfem::Operator {
public:
  HertzOperator(mfem::ParMesh &pmesh, int order,
                mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
                mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                hephaestus::BCMap &bc_map,
                hephaestus::DomainProperties &domain_properties,
                hephaestus::Sources &sources,
                hephaestus::InputParameters &solver_options)
      : myid_(0), num_procs_(1), _order(order), pmesh_(&pmesh),
        _fespaces(fespaces), _variables(variables), _bc_map(bc_map),
        _sources(sources), _domain_properties(domain_properties),
        _solver_options(solver_options){};

  ~HertzOperator(){};

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
  mfem::ParComplexGridFunction *e_;
  mfem::ParComplexLinearForm *jd_;
  mfem::ParSesquilinearForm *a1_;
  mfem::ComplexOperator::Convention conv_ = mfem::ComplexOperator::HERMITIAN;

  int myid_;
  int num_procs_;
  const int _order;
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

  mfem::Coefficient *muInvCoef_; // Dia/Paramagnetic Material Coefficient
  mfem::Coefficient *massCoef_;  // -omega^2 epsilon
  mfem::Coefficient *lossCoef_;  // omega sigma

  mfem::VectorCoefficient *jrCoef_; // Volume Current Density Function
  mfem::VectorCoefficient *jiCoef_; // Volume Current Density Function

  mfem::Array<int> ess_bdr_tdofs_;
};

//
// Specifies output interfaces of a time-domain EM formulation.
class HertzFormulation {
  // std::vector<mfem::ParGridFunction *> local_trial_vars, local_test_vars;
  std::string frequency_coef_name, permittivity_coef_name,
      reluctivity_coef_name, conductivity_coef_name, h_curl_var_name;

public:
  hephaestus::HertzOperator *fd_operator;
  mfem::ConstantCoefficient oneCoef;
  mfem::ConstantCoefficient *freqCoef;
  HertzFormulation();

  // virtual hephaestus::TimeDependentEquationSystem *CreateEquationSystem();

  virtual hephaestus::HertzOperator *CreateFrequencyDomainOperator(
      mfem::ParMesh &pmesh, int order,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options);

  virtual void RegisterMissingVariables(
      mfem::ParMesh &pmesh,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables);

  virtual void
  RegisterAuxKernels(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                     hephaestus::AuxKernels &auxkernels){};

  virtual void
  RegisterCoefficients(hephaestus::DomainProperties &domain_properties);
};
} // namespace hephaestus
