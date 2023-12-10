#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "problem_builder_base.hpp"

namespace hephaestus {
class ProblemOperator : public mfem::Operator {
public:
  ProblemOperator(hephaestus::Problem &problem) : _problem(problem){};

  ~ProblemOperator(){};

  virtual void SetGridFunctions();
  virtual void Init(mfem::Vector &X);
  virtual void Solve(mfem::Vector &X){};
  virtual mfem::Solver *getJacobianSolver() {
    return _problem._jacobian_solver.get();
  };
  virtual mfem::NewtonSolver *getNonlinearSolver() {
    return _problem._nonlinear_solver.get();
  };
  void Mult(const mfem::Vector &x, mfem::Vector &y) const override{};

  mfem::Array<int> true_offsets, block_trueOffsets;

  // Vector of names of state gridfunctions used in formulation, ordered by
  // appearance in block vector during solve.
  std::vector<std::string> trial_var_names;
  std::vector<mfem::ParGridFunction *> trial_variables;

  mfem::BlockVector trueX, trueRhs;
  mfem::OperatorHandle _equation_system_operator;

protected:
  hephaestus::Problem &_problem;
};

} // namespace hephaestus
