#pragma once
#include "vector_coefficient_aux.hpp"

// Specify postprocessors that depend on one or more gridfunctions
namespace hephaestus {

void cross_product(mfem::Vector &va, mfem::Vector &vb, mfem::Vector &V);

// The LorentzForceDensityAuxCoefficient object will contain a reference to the
// current density and magnetic flux density grid functions,
// and returns the Lorentz force density coefficient.
class LorentzForceDensityAuxCoefficient : public mfem::VectorCoefficient {
private:
  mfem::ParGridFunction &J_gf;
  mfem::ParGridFunction &B_gf;

public:
  LorentzForceDensityAuxCoefficient(mfem::ParGridFunction &J_gf_,
                                    mfem::ParGridFunction &B_gf_)
      : mfem::VectorCoefficient(3), J_gf(J_gf_), B_gf(B_gf_) {}
  virtual void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
                    const mfem::IntegrationPoint &ip);
  virtual ~LorentzForceDensityAuxCoefficient() {}
};

// Project a stored scalar Coefficient onto a (scalar) GridFunction
class LorentzForceDensityAux : public VectorCoefficientAux {
private:
  mfem::ParGridFunction *j_gf;
  mfem::ParGridFunction *b_gf;

  const std::string _magnetic_flux_density_name;
  const std::string _current_density_name;

public:
  LorentzForceDensityAux(const std::string &lorentz_force_density_gf_name,
                         const std::string &lorentz_force_density_coef_name,
                         const std::string &magnetic_flux_density_name,
                         const std::string &current_density_name);

  void Init(const hephaestus::GridFunctions &gridfunctions,
            hephaestus::Coefficients &coefficients) override;
};

} // namespace hephaestus
