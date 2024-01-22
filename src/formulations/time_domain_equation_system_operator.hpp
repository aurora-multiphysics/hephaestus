#pragma once
#include "../common/pfem_extras.hpp"
#include "auxsolvers.hpp"
#include "equation_system.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus
{

std::string GetTimeDerivativeName(const std::string & name);

std::vector<std::string> GetTimeDerivativeNames(std::vector<std::string> gridfunction_names);

// Specifies output interfaces of a time-domain EM formulation.
class TimeDomainEquationSystemOperator : public mfem::TimeDependentOperator
{
public:
  TimeDomainEquationSystemOperator(mfem::ParMesh & pmesh,
                                   hephaestus::FESpaces & fespaces,
                                   hephaestus::GridFunctions & gridfunctions,
                                   hephaestus::BCMap & bc_map,
                                   hephaestus::Coefficients & coefficients,
                                   hephaestus::Sources & sources,
                                   hephaestus::InputParameters & solver_options)
    : pmesh_(&pmesh),
      _fespaces(fespaces),
      _gridfunctions(gridfunctions),
      _bc_map(bc_map),
      _sources(sources),
      _coefficients(coefficients),
      _solver_options(solver_options){};

  ~TimeDomainEquationSystemOperator() override = default;

  virtual void SetGridFunctions();
  virtual void Init(mfem::Vector & X);
  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override;
  void SetEquationSystem(hephaestus::TimeDependentEquationSystem * equation_system);

  mfem::Array<int> true_offsets, block_trueOffsets;
  // Vector of names of state gridfunctions used in formulation, ordered by
  // appearance in block vector during solve.
  std::vector<std::string> state_var_names;
  // Vector of names of recognised auxiliary gridfunctions that can be
  // calculated from formulation,
  std::vector<std::string> aux_var_names;
  // Vector of names of active auxiliary gridfunctions that are being calculated
  // in formulation,
  std::vector<std::string> active_aux_var_names;

  std::vector<mfem::ParGridFunction *> local_trial_vars, local_test_vars;

  hephaestus::TimeDependentEquationSystem * _equation_system;

  int myid_{0};
  int num_procs_{1};
  mfem::ParMesh * pmesh_;
  hephaestus::FESpaces & _fespaces;
  hephaestus::GridFunctions & _gridfunctions;
  hephaestus::BCMap & _bc_map;
  hephaestus::Sources & _sources;
  hephaestus::Coefficients & _coefficients;
  hephaestus::InputParameters & _solver_options;

  mutable std::unique_ptr<hephaestus::DefaultGMRESSolver> solver{nullptr};
  mutable std::unique_ptr<hephaestus::DefaultHCurlPCGSolver> jacobian_solver{nullptr};

  mfem::OperatorHandle blockA;
  mfem::BlockVector trueX, trueRhs;
};

} // namespace hephaestus
