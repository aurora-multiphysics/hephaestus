#pragma once
#include "coefficients.hpp"
#include "gridfunctions.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "mfem.hpp"

// Specify classes that perform auxiliary calculations on GridFunctions or
// Coefficients.
namespace hephaestus {

class AuxSolver {
public:
  AuxSolver() = default;
  int priority{0};

  virtual void Init(const hephaestus::GridFunctions &gridfunctions,
                    hephaestus::Coefficients &coefficients) = 0;

  virtual void Solve(double t = 0.0) = 0;
  // Set priority. Lower values are evaluated first.
  void SetPriority(int prty) { priority = prty; };
};

struct AuxCompare {
  bool const operator()(AuxSolver *lhs, AuxSolver *rhs) const {
    return (lhs->priority) < (rhs->priority);
  }
};

} // namespace hephaestus
