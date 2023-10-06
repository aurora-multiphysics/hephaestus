#pragma once
#include "../common/pfem_extras.hpp"
#include "auxsolvers.hpp"
#include "equation_system.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus {
class SteadyStateEquationSystemOperator : public mfem::Operator {
public:
  SteadyStateEquationSystemOperator(
      mfem::ParMesh &pmesh, hephaestus::FESpaces &fespaces,
      hephaestus::GridFunctions &gridfunctions, hephaestus::BCMap &bc_map,
      hephaestus::Coefficients &coefficients, hephaestus::Sources &sources,
      hephaestus::InputParameters &solver_options)
      : myid_(0), num_procs_(1), pmesh_(&pmesh), _fespaces(fespaces),
        _gridfunctions(gridfunctions), _bc_map(bc_map), _sources(sources),
        _coefficients(coefficients), _solver_options(solver_options){};

  ~SteadyStateEquationSystemOperator(){};

  virtual void SetGridFunctions();
  virtual void Init(mfem::Vector &X);
  virtual void Solve(mfem::Vector &X);
  void Mult(const mfem::Vector &x, mfem::Vector &y) const override{};

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

  int myid_;
  int num_procs_;
  mfem::ParMesh *pmesh_;
  hephaestus::FESpaces &_fespaces;
  hephaestus::GridFunctions &_gridfunctions;
  hephaestus::BCMap &_bc_map;
  hephaestus::Sources &_sources;
  hephaestus::Coefficients &_coefficients;
  hephaestus::InputParameters &_solver_options;

  mfem::OperatorHandle blockA;
  mfem::BlockVector trueX, trueRhs;
};

} // namespace hephaestus