#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "problem_builder_base.hpp"

namespace hephaestus {

std::string GetTimeDerivativeName(const std::string &name);

std::vector<std::string>
GetTimeDerivativeNames(std::vector<std::string> gridfunction_names);

// Specifies output interfaces of a time-domain EM formulation.
class TimeDomainProblemOperator : public mfem::TimeDependentOperator {
public:
  TimeDomainProblemOperator(hephaestus::Problem &problem) : _problem(problem) {
    ;
  };

  ~TimeDomainProblemOperator(){};

  virtual void SetGridFunctions();
  virtual void Init(mfem::Vector &X);
  virtual void ImplicitSolve(const double dt, const mfem::Vector &X,
                             mfem::Vector &dX_dt) override;

  virtual void buildEquationSystemOperator(double dt);

  void
  SetEquationSystem(hephaestus::TimeDependentEquationSystem *equation_system);

  virtual mfem::Solver *getJacobianSolver() {
    return _problem._jacobian_solver.get();
  };
  virtual mfem::NewtonSolver *getNonlinearSolver() {
    return _problem._nonlinear_solver.get();
  };

  mfem::Array<int> true_offsets, block_trueOffsets;

  // Vector of names of state gridfunctions used in formulation, ordered by
  // appearance in block vector during solve.
  std::vector<std::string> trial_var_names;
  std::vector<mfem::ParGridFunction *> trial_variable_time_derivatives,
      trial_variables;

  hephaestus::TimeDependentEquationSystem *_equation_system;

  mfem::BlockVector trueX, trueRhs;
  mfem::OperatorHandle _equation_system_operator;

protected:
  hephaestus::Problem &_problem;
};

} // namespace hephaestus
