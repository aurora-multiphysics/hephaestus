#pragma once
#include "problem_builder_base.hpp"
#include "time_domain_problem_operator.hpp"

namespace hephaestus
{
// Stores data required to describe a time domain formulation
using TimeDomainEquationSystemProblem =
    EquationSystemProblemTemplate<TimeDomainEquationSystemProblemOperator>;

// Builder class of a time-domain EM formulation.
class TimeDomainEquationSystemProblemBuilder
  : public EquationSystemProblemBuilder<TimeDomainEquationSystemProblem>
{
public:
  TimeDomainEquationSystemProblemBuilder()
    : _problem(std::make_unique<hephaestus::TimeDomainEquationSystemProblem>())
  {
  }

  ~TimeDomainEquationSystemProblemBuilder() override = default;

  virtual std::unique_ptr<hephaestus::TimeDomainEquationSystemProblem> ReturnProblem()
  {
    return std::move(_problem);
  }

  static std::vector<mfem::ParGridFunction *>
  RegisterTimeDerivatives(std::vector<std::string> gridfunction_names,
                          hephaestus::GridFunctions & gridfunctions);

  void RegisterFESpaces() override {}

  void RegisterGridFunctions() override;

  void RegisterAuxSolvers() override {}

  void RegisterCoefficients() override {}

  void SetOperatorGridFunctions() override;

  void ConstructOperator() override;

  void ConstructState() override;

  void ConstructTimestepper() override;

  void InitializeKernels() override
  {
    ProblemBuilder::InitializeKernels();

    GetProblem()->GetEquationSystem()->Init(GetProblem()->_gridfunctions,
                                            GetProblem()->_fespaces,
                                            GetProblem()->_bc_map,
                                            GetProblem()->_coefficients);
  }

protected:
  std::unique_ptr<hephaestus::TimeDomainEquationSystemProblem> _problem{nullptr};
  mfem::ConstantCoefficient _one_coef{1.0};

  TimeDomainEquationSystemProblem * GetProblem() override { return _problem.get(); };
};

} // namespace hephaestus
