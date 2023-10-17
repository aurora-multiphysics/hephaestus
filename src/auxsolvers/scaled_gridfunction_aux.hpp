#pragma once
#include "auxsolvers.hpp"

namespace hephaestus {

class ScaledGridFunctionAuxSolver : public AuxSolver {
private:
  // ConstantCoefficient to scale input gridfunction by
  mfem::ConstantCoefficient *coef;
  // Input gridfunction to be scaled by a scalar coefficient
  mfem::ParGridFunction *input_gf;

  // Gridfunction in which to store result
  mfem::ParGridFunction *scaled_gf;

  std::string coef_name;
  std::string input_gf_name;
  std::string scaled_gf_name;

public:
  ScaledGridFunctionAuxSolver(const hephaestus::InputParameters &params);

  void Init(const hephaestus::GridFunctions &gridfunctions,
            hephaestus::Coefficients &coefficients) override;

  void Solve(double t = 0.0) override;
  void buildHCurlMass();
  void buildWeakDiv();
};

} // namespace hephaestus
