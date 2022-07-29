#pragma once
#include "materials.hpp"
#include "mfem.hpp"
#include "variables.hpp"

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

class Postprocessor {
public:
  Postprocessor() {}

  Postprocessor(const hephaestus::InputParameters &params){};

  virtual void
  Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       hephaestus::DomainProperties &domain_properties){};

  virtual void Update(double t = 0.0){};
};

class Postprocessors : public mfem::NamedFieldsMap<hephaestus::Postprocessor> {
public:
  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::DomainProperties &domain_properties);
  void Update(double t = 0.0);
};

// Class to calculate and store the L2 error
// of a grid function with respect to a (Vector)Coefficient
class L2ErrorVectorPostprocessor : public Postprocessor {

public:
  L2ErrorVectorPostprocessor(){};
  L2ErrorVectorPostprocessor(const hephaestus::InputParameters &params);

  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::DomainProperties &domain_properties) override;

  virtual void Update(double t = 0.0) override;

  std::string var_name;      // name of the variable
  std::string vec_coef_name; // name of the vector coefficient

  mfem::Array<double> times;
  mfem::Array<HYPRE_BigInt> ndofs;
  mfem::Array<double> l2_errs;
  mfem::ParGridFunction *gf;
  mfem::VectorCoefficient *vec_coeff;
};

} // namespace hephaestus
