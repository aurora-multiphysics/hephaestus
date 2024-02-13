#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "problem_builder_base.hpp"

namespace hephaestus
{

class ProblemOperator : public mfem::Operator
{
public:
  ProblemOperator(hephaestus::Problem & problem) : _problem(problem){};

  virtual void SetGridFunctions();
  virtual void Init(mfem::Vector & X);
  virtual void Solve(mfem::Vector & X){};
  void Mult(const mfem::Vector & x, mfem::Vector & y) const override {}

  mfem::Array<int> _true_offsets, _block_true_offsets;

  // Vector of names of state gridfunctions used in formulation, ordered by
  // appearance in block vector during solve.
  std::vector<std::string> _trial_var_names;
  std::vector<mfem::ParGridFunction *> _trial_variables;

  mfem::BlockVector _true_x, _true_rhs;
  mfem::OperatorHandle _equation_system_operator;

  void SetEquationSystem(std::unique_ptr<EquationSystem> new_equation_system)
  {
    _equation_system.reset();
    _equation_system = std::move(new_equation_system);
  }

  inline EquationSystem * GetEquationSystem() const
  {
    if (!_equation_system)
    {
      MFEM_ABORT("No equation system has been added to TimeDomainProblemOperator.");
    }

    return _equation_system.get();
  }

protected:
  hephaestus::Problem & _problem;
  std::unique_ptr<EquationSystem> _equation_system{nullptr};
};

} // namespace hephaestus