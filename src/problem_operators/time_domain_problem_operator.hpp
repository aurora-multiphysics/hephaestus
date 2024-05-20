#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "problem_builder_base.hpp"
#include "problem_operator_base.hpp"

namespace hephaestus
{

std::string GetTimeDerivativeName(const std::string & name);

std::vector<std::string> GetTimeDerivativeNames(std::vector<std::string> gridfunction_names);

/// Problem operator for time-dependent problems with no equation system. The user will need to subclass this since the solve is not
/// implemented.
class TimeDomainProblemOperator : public mfem::TimeDependentOperator, public ProblemOperatorBase
{
public:
  TimeDomainProblemOperator(hephaestus::Problem & problem) : ProblemOperatorBase(problem) {}
  ~TimeDomainProblemOperator() override = default;

  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override {}

protected:
  int & Width() final { return mfem::TimeDependentOperator::width; }
  int & Height() final { return mfem::TimeDependentOperator::height; }
};

} // namespace hephaestus