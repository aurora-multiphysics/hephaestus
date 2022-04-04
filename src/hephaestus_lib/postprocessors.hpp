#pragma once
#include "mfem.hpp"

namespace hephaestus {

void cross_product(mfem::Vector &va, mfem::Vector &vb, mfem::Vector &V);

// The JouleHeatingCoefficient object will contain a reference to the electric
// field grid function, and the conductivity sigma, and returns sigma E dot E at
// a point.
class JouleHeatingCoefficient : public mfem::Coefficient {
private:
  mfem::ParGridFunction &E_gf;
  mfem::PWCoefficient sigma;

public:
  JouleHeatingCoefficient(const mfem::PWCoefficient &sigma_,
                          mfem::ParGridFunction &E_gf_)
      : E_gf(E_gf_), sigma(sigma_) {}
  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip);
  virtual ~JouleHeatingCoefficient() {}
};

// The LorentzForceVectorCoefficient object will contain a reference to the
// electric field grid function, and the conductivity sigma, and returns sigma E
// dot E at a point.
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
