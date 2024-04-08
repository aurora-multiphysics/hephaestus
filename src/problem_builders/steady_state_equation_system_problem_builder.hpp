#pragma once
#include "problem_operator.hpp"
#include "problem_builder_base.hpp"
namespace hephaestus
{

/// Class for steady-state problems with an equation system.
using SteadyStateEquationSystemProblem =
    EquationSystemProblemTemplate<EquationSystemProblemOperator>;

// Builder class of a frequency-domain problem.
class SteadyStateEquationSystemProblemBuilder
  : public EquationSystemProblemBuilder<SteadyStateEquationSystemProblem>
{
public:
  SteadyStateEquationSystemProblemBuilder()
    : _problem(std::make_unique<hephaestus::SteadyStateEquationSystemProblem>())
  {
  }

  ~SteadyStateEquationSystemProblemBuilder() override = default;

  virtual std::unique_ptr<hephaestus::SteadyStateEquationSystemProblem> ReturnProblem()
  {
    return std::move(_problem);
  }

  void RegisterFESpaces() override {}

  void RegisterGridFunctions() override {}

  void RegisterAuxSolvers() override {}

  void RegisterCoefficients() override {}

  void SetOperatorGridFunctions() override;

  void ConstructOperator() override;

  void ConstructState() override;

  void ConstructTimestepper() override {}

  void InitializeKernels() override
  {
    ProblemBuilder::InitializeKernels();

    GetProblem()->GetEquationSystem()->Init(GetProblem()->_gridfunctions,
                                            GetProblem()->_fespaces,
                                            GetProblem()->_bc_map,
                                            GetProblem()->_coefficients);
  }

protected:
  std::unique_ptr<hephaestus::SteadyStateEquationSystemProblem> _problem{nullptr};
  mfem::ConstantCoefficient _one_coef{1.0};

  SteadyStateEquationSystemProblem * GetProblem() override { return _problem.get(); };
};

} // namespace hephaestus
