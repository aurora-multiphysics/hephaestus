#pragma once
#include "auxsolver_base.hpp"

// Specify classes that perform auxiliary calculations on GridFunctions or
// Coefficients.
namespace hephaestus {

// Class to allow creation of scalar coefficients that are coupled to
// gridfunctions. Defines a Coefficient that is updated from a GridFunction (to
// be moved into Coefficients?)
class CoupledCoefficient : public mfem::Coefficient, public AuxSolver {
protected:
  // pointer to coupled variable (could just be to vars?)
  mfem::ParGridFunction *gf;
  double scalar_val;

public:
  CoupledCoefficient(const hephaestus::InputParameters &params);

  virtual void
  Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       hephaestus::Coefficients &domain_properties) override;

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip) override;

  void Solve(double t = 0.0) override{};
  virtual ~CoupledCoefficient() {}

  std::string coupled_var_name; // name of the variable
};

} // namespace hephaestus
