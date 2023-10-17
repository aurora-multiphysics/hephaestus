#pragma once
#include "auxsolver_base.hpp"

namespace hephaestus {

// Scale a gridfunction in H(Curl) or H(Div) by a scalar Coefficient, and store
// the result. Suitable for solving for H(Div) or H(Curl) conforming fields for
// expressions like J = σE
class ScaledVectorGridFunctionAux : public AuxSolver {
public:
  ScaledVectorGridFunctionAux(const hephaestus::InputParameters &params);
  virtual ~ScaledVectorGridFunctionAux();
  virtual void Init(const hephaestus::GridFunctions &gridfunctions,
                    hephaestus::Coefficients &coefficients) override;
  virtual void buildBilinearForm();
  virtual void buildMixedBilinearForm();
  virtual void Solve(double t = 0.0) override;

protected:
  // Pointers to store trial and test FE spaces
  mfem::ParFiniteElementSpace *trial_fes;
  mfem::ParFiniteElementSpace *test_fes;

  // Bilinear forms
  mfem::ParBilinearForm *a;
  mfem::ParMixedBilinearForm *a_mixed;

  // Coefficient to scale input gridfunction by
  mfem::Coefficient *coef;

private:
  std::string coef_name;
  std::string input_gf_name;
  std::string scaled_gf_name;
  const hephaestus::InputParameters solver_options;

  // Input gridfunction to be scaled by a scalar coefficient
  mfem::ParGridFunction *input_gf;

  // Gridfunction in which to store result
  mfem::ParGridFunction *scaled_gf;

  // Operator matrices
  mfem::HypreParMatrix *a_mat;
  mfem::HypreParMatrix *mixed_mat;

  // Solver
  hephaestus::DefaultH1PCGSolver *solver;
};
} // namespace hephaestus
