#pragma once
#include "problem_operator.hpp"
#include "problem_builder_base.hpp"
namespace hephaestus
{

/// Class for steady-state problems with no equation system.
using SteadyStateProblem = ProblemTemplate<ProblemOperator>;

class SteadyStateProblemBuilder : public hephaestus::ProblemBuilder
{
public:
  SteadyStateProblemBuilder() : _problem(std::make_unique<hephaestus::SteadyStateProblem>()) {}

  ~SteadyStateProblemBuilder() override = default;

  virtual std::unique_ptr<hephaestus::SteadyStateProblem> ReturnProblem()
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

protected:
  std::unique_ptr<hephaestus::SteadyStateProblem> _problem{nullptr};
  mfem::ConstantCoefficient _one_coef{1.0};

  hephaestus::SteadyStateProblem * GetProblem() override { return _problem.get(); };
};

} // namespace hephaestus
