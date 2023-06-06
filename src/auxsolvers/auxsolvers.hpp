#pragma once
#include "inputs.hpp"
#include "materials.hpp"
#include "mfem.hpp"

// Specify kernels to modify and solve auxiliary variables or coefficients using
// one or more DomainProperties.
namespace hephaestus {

class AuxSolver {
public:
  AuxSolver() = default;

  virtual void
  Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       hephaestus::DomainProperties &domain_properties) = 0;

  virtual void Solve(double t = 0.0) = 0;
};

class AuxSolvers : public mfem::NamedFieldsMap<hephaestus::AuxSolver> {
public:
  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::DomainProperties &domain_properties);
  void Solve(double t = 0.0);
};

class CurlAuxSolver : public AuxSolver {
public:
  CurlAuxSolver(const hephaestus::InputParameters &params);

  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::DomainProperties &domain_properties) override;

  void Solve(double t = 0.0) override;

  std::string var_name;      // name of the variable
  std::string curl_var_name; // Variable in which to store curl

  mfem::ParGridFunction *u_, *curl_u_;
  mfem::ParDiscreteLinearOperator *curl;
};

class CoefficientAuxSolver : public AuxSolver {
public:
  CoefficientAuxSolver(const hephaestus::InputParameters &params);

  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::DomainProperties &domain_properties) override;

  void Solve(double t = 0.0) override;

  std::string var_name;  // name of the variable
  std::string coef_name; // name of the coefficient

  mfem::ParGridFunction *gf;
  mfem::Coefficient *coeff;
};

class VectorCoefficientAuxSolver : public AuxSolver {
public:
  VectorCoefficientAuxSolver(const hephaestus::InputParameters &params);

  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::DomainProperties &domain_properties) override;

  void Solve(double t = 0.0) override;

  std::string var_name;      // name of the variable
  std::string vec_coef_name; // name of the vector coefficient

  mfem::ParGridFunction *gf;
  mfem::VectorCoefficient *vec_coeff;
};

// Class to allow creation of coefficients that are coupled to gridfunctions.
// Should be separate to AuxSolvers? Most similar to MOOSE materials...
class CoupledCoefficient : public mfem::Coefficient, public AuxSolver {
protected:
  // pointer to coupled variable (could just be to vars?)
  mfem::ParGridFunction *gf;
  double scalar_val;

public:
  CoupledCoefficient(const hephaestus::InputParameters &params);

  virtual void
  Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       hephaestus::DomainProperties &domain_properties) override;

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip) override;

  void Solve(double t = 0.0) override{};
  virtual ~CoupledCoefficient() {}

  std::string coupled_var_name; // name of the variable
};

} // namespace hephaestus
