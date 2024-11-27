#pragma once
#include "auxsolver_base.hpp"

// Specify postprocessors that depend on one or more gridfunctions
namespace hephaestus
{

double calcSurfaceForceDensity(mfem::ParGridFunction * b_field,
                               mfem::ParGridFunction * h_field,
                               int attr,
                               mfem::Coefficient & q,
                               mfem::Coefficient & mu,
                               bool use_face_attr);

// Class to calculate and store the flux of a vector GridFunction through a surface
// at each timestep, optionally scaled by a coefficient.
class MagneticForceAux : public AuxSolver
{

public:
  MagneticForceAux() = default;
  MagneticForceAux(std::string b_name,
                   std::string h_name,
                   mfem::Array<int> attr,
                   std::string coef_name = "",
                   bool use_face_attr = true);

  ~MagneticForceAux() override = default;

  void Init(const hephaestus::GridFunctions & gridfunctions,
            hephaestus::Coefficients & coefficients) override;

  void Solve(double t = 0.0) override;

  // void WriteForces(std::string fname, mfem::ParGridFunction & gf, int attr);

  std::string _b_name;    // name of the vector variable
  std::string _h_name;    // name of the vector variable
  std::string _coef_name; // name of the coefficient

  mfem::Array<double> _times;
  mfem::Array<int> _attr_out;
  mfem::Array<double> _forces;

  mfem::ParGridFunction * _b_gf{nullptr};
  mfem::ParGridFunction * _h_gf{nullptr};
  mfem::ParGridFunction * _gf{nullptr};
  mfem::Coefficient * _mu_coef{nullptr};

  mfem::Array<int> _attr;
  bool _use_face_attr;
};

} // namespace hephaestus
