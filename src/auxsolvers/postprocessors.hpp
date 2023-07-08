#pragma once
#include "auxsolvers.hpp"

// Specify postprocessors that depend on one or more gridfunctions
namespace hephaestus {

void cross_product(mfem::Vector &va, mfem::Vector &vb, mfem::Vector &V);

// The LorentzForceVectorCoefficient object will contain a reference to the
// current density and magnetic flux density grid functions,
// and returns the Lorentz force density coefficient.
class LorentzForceVectorCoefficient : public mfem::VectorCoefficient {
private:
  mfem::ParGridFunction &J_gf;
  mfem::ParGridFunction &B_gf;

public:
  LorentzForceVectorCoefficient(mfem::ParGridFunction &J_gf_,
                                mfem::ParGridFunction &B_gf_)
      : mfem::VectorCoefficient(3), J_gf(J_gf_), B_gf(B_gf_) {}
  virtual void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
                    const mfem::IntegrationPoint &ip);
  virtual ~LorentzForceVectorCoefficient() {}
};

} // namespace hephaestus
