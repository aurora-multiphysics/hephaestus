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
  TimeDomainEquationSystemProblemBuilder() = default;
  ~TimeDomainEquationSystemProblemBuilder() override = default;

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
  mfem::ConstantCoefficient _one_coef{1.0};
};

} // namespace hephaestus
