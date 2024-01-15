#pragma once
#include "coefficients.hpp"
#include "gridfunctions.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "mfem.hpp"

// Specify classes that perform auxiliary calculations on GridFunctions or
// Coefficients.
namespace hephaestus
{

class AuxSolver;

class AuxSolver
{
public:
  AuxSolver() = default;

  virtual void Init(const hephaestus::GridFunctions & gridfunctions,
                    hephaestus::Coefficients & coefficients) = 0;

  virtual void Solve(double t = 0.0) = 0;

  // Set priority. Lower values are evaluated first.
  void SetPriority(const int priority) { _priority = priority; };

  inline int priority() const { return _priority; }

  // Priority comparator.
  static bool comparator(const AuxSolver * first, const AuxSolver * second)
  {
    return (first->priority() < second->priority());
  }

private:
  int _priority{0};
};

} // namespace hephaestus
