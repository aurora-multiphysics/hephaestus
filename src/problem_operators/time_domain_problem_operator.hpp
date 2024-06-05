#pragma once
#include "../common/pfem_extras.hpp"
#include "hephaestus_solvers.hpp"
#include "problem_operator_base.hpp"

namespace hephaestus
{

std::string GetTimeDerivativeName(const std::string & name);

std::vector<std::string> GetTimeDerivativeNames(std::vector<std::string> gridfunction_names);

/// Problem operator for time-dependent problems with no equation system. The user will need to subclass this since the solve is not
/// implemented.
class TimeDomainProblemOperator : virtual public mfem::TimeDependentOperator, public ProblemOperator
{
public:
  TimeDomainProblemOperator(hephaestus::Problem & problem);
  ~TimeDomainProblemOperator() override = default;

  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override {}

  /// Wrapper around the ODE solver's Step method using the block vector.
  virtual void Step(mfem::real_t & t, mfem::real_t & dt)
  {
    _ode_solver->Step(*_block_vector, t, dt);
  }

  void Init() override;

  void Update() override;

protected:
  int & Width() final { return mfem::TimeDependentOperator::width; }
  int & Height() final { return mfem::TimeDependentOperator::height; }

  /// Construct the ODE solver. Called in constructor. Override in derived classes.
  virtual void ConstructTimestepper();

private:
  /// Store the ODE solver.
  std::unique_ptr<mfem::ODESolver> _ode_solver{nullptr};
};

} // namespace hephaestus