#pragma once
#include "auxsolver_base.hpp"

// Specify postprocessors that depend on one or more variables
namespace hephaestus {

// Class to calculate and store the L2 error
// of a grid function with respect to a (Vector)Coefficient
class L2ErrorVectorPostprocessor : public AuxSolver {

public:
  L2ErrorVectorPostprocessor(){};
  L2ErrorVectorPostprocessor(const hephaestus::InputParameters &params);

  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::Coefficients &coefficients) override;

  virtual void Solve(double t = 0.0) override;

  std::string var_name;      // name of the variable
  std::string vec_coef_name; // name of the vector coefficient

  mfem::Array<double> times;
  mfem::Array<HYPRE_BigInt> ndofs;
  mfem::Array<double> l2_errs;
  mfem::ParGridFunction *gf;
  mfem::VectorCoefficient *vec_coeff;
};

} // namespace hephaestus
