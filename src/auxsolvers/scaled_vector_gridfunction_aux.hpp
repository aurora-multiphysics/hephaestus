#pragma once
#include "auxsolver_base.hpp"

namespace hephaestus {

// Scale a gridfunction in H(Curl) or H(Div) by a scalar Coefficient, and store
// the result. Suitable for solving for H(Div) or H(Curl) conforming fields for
// expressions like v = a*σ*u
class ScaledVectorGridFunctionAux : public AuxSolver {
public:
  ScaledVectorGridFunctionAux(
      const std::string &input_gf_name, const std::string &scaled_gf_name,
      const std::string &coef_name, const double &aConst = 1.0,
      const hephaestus::InputParameters &solver_options =
          hephaestus::InputParameters());
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
  // Optional constant to scale input gridfunction by

private:
  const std::string _input_gf_name;
  const std::string _scaled_gf_name;
  const std::string _coef_name;
  const double _aConst;
  const hephaestus::InputParameters _solver_options;

  // Input gridfunction to be scaled by a scalar coefficient
  mfem::ParGridFunction *input_gf;

  // Gridfunction in which to store result
  mfem::ParGridFunction *scaled_gf;

  // Operator matrices
  mfem::HypreParMatrix *a_mat;
  mfem::HypreParMatrix *mixed_mat;

  // Solver
  hephaestus::DefaultJacobiPCGSolver *solver;
};
} // namespace hephaestus
