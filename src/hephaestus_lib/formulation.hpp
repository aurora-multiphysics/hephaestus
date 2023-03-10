#pragma once
#include "../common/pfem_extras.hpp"
#include "auxkernels.hpp"
#include "equation_system.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {

template <typename T>
std::vector<T *>
populateVectorFromNamedFieldsMap(mfem::NamedFieldsMap<T> nfmap,
                                 std::vector<std::string> keys) {
  std::vector<T *> result;
  for (auto &key : keys) {
    if (nfmap.Has(key)) {
      result.push_back(nfmap.Get(key));
    } else {
      MFEM_ABORT("Key " << key << " not found in NamedFieldsMap");
    }
  }
  return result;
};

std::string GetTimeDerivativeName(std::string name);

std::vector<std::string>
GetTimeDerivativeNames(std::vector<std::string> gridfunction_names);

// Specifies output interfaces of a time-domain EM formulation.
class TimeDomainEquationSystemOperator : public mfem::TimeDependentOperator {
public:
  TimeDomainEquationSystemOperator(
      mfem::ParMesh &pmesh, int order,
      mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
      hephaestus::BCMap &bc_map,
      hephaestus::DomainProperties &domain_properties,
      hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
      : myid_(0), num_procs_(1), _order(order), pmesh_(&pmesh),
        _fespaces(fespaces), _variables(variables), _bc_map(bc_map),
        _sources(sources), _domain_properties(domain_properties),
        _solver_options(solver_options){};

  ~TimeDomainEquationSystemOperator(){};

  virtual void Init(mfem::Vector &X);
  virtual void ImplicitSolve(const double dt, const mfem::Vector &X,
                             mfem::Vector &dX_dt) override;
  void
  SetEquationSystem(hephaestus::TimeDependentEquationSystem *equation_system);

  virtual void SetVariables();

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
};

//
// Specifies output interfaces of a time-domain EM formulation.
class TransientFormulation {
  // std::vector<mfem::ParGridFunction *> local_trial_vars, local_test_vars;

public:
  hephaestus::TimeDependentEquationSystem *equation_system;
  hephaestus::TimeDomainEquationSystemOperator *td_operator;
  mfem::ConstantCoefficient oneCoef;

  TransientFormulation();

  virtual hephaestus::TimeDependentEquationSystem *CreateEquationSystem();

  virtual hephaestus::TimeDomainEquationSystemOperator *
  CreateTimeDomainOperator(
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
      mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) {
    if (equation_system != NULL) {
      TransientFormulation::RegisterTimeDerivatives(equation_system->var_names,
                                                    variables);
    }
  };

  static std::vector<mfem::ParGridFunction *> RegisterTimeDerivatives(
      std::vector<std::string> gridfunction_names,
      mfem::NamedFieldsMap<mfem::ParGridFunction> &gridfunctions);

  virtual void
  RegisterAuxKernels(mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                     hephaestus::AuxKernels &auxkernels){};

  virtual void
  RegisterCoefficients(hephaestus::DomainProperties &domain_properties){};
};
} // namespace hephaestus
