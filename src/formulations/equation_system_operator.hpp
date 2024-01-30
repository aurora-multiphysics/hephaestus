#pragma once
#include "../common/pfem_extras.hpp"
#include "auxsolvers.hpp"
#include "equation_system.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "sources.hpp"

namespace hephaestus
{
class EquationSystemOperator : public mfem::Operator
{
public:
  EquationSystemOperator(mfem::ParMesh & pmesh,
                         hephaestus::FESpaces & fespaces,
                         hephaestus::GridFunctions & gridfunctions,
                         hephaestus::BCMap & bc_map,
                         hephaestus::Coefficients & coefficients,
                         hephaestus::Sources & sources,
                         hephaestus::InputParameters & solver_options)
    : _pmesh(&pmesh),
      _fespaces(fespaces),
      _gridfunctions(gridfunctions),
      _bc_map(bc_map),
      _sources(sources),
      _coefficients(coefficients),
      _solver_options(solver_options)
  {
  }

  ~EquationSystemOperator() override = default;

  virtual void SetGridFunctions();
  virtual void Init(mfem::Vector & X);
  virtual void Solve(mfem::Vector & X);
  void Mult(const mfem::Vector & x, mfem::Vector & y) const override {}
  void SetSolver(const mfem::HypreParMatrix & M);
  void SetSolver(const mfem::HypreParMatrix & M, mfem::ParFiniteElementSpace * edge_fespace);

  mfem::Array<int> _true_offsets, _block_true_offsets;
  // Vector of names of state gridfunctions used in formulation, ordered by
  // appearance in block vector during solve.
  std::vector<std::string> _state_var_names;
  // Vector of names of recognised auxiliary gridfunctions that can be
  // calculated from formulation,
  std::vector<std::string> _aux_var_names;
  // Vector of names of active auxiliary gridfunctions that are being calculated
  // in formulation,
  std::vector<std::string> _active_aux_var_names;

  std::vector<mfem::ParGridFunction *> _local_test_vars;

  int _myid{0};
  int _num_procs{1};
  mfem::ParMesh * _pmesh;
  hephaestus::FESpaces & _fespaces;
  hephaestus::GridFunctions & _gridfunctions;
  hephaestus::BCMap & _bc_map;
  hephaestus::Sources & _sources;
  hephaestus::Coefficients & _coefficients;
  hephaestus::InputParameters & _solver_options;

  mfem::OperatorHandle _block_a;
  mfem::BlockVector _true_x, _true_rhs;

  std::unique_ptr<mfem::HypreSolver> _solver{nullptr};
};

} // namespace hephaestus
