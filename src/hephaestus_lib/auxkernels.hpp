#pragma once
#include "materials.hpp"
#include "mfem.hpp"
#include "variables.hpp"

// Specify kernels to modify and solve auxiliary variables using one or more
// DomainProperties.
namespace hephaestus {

class AuxKernel {
public:
  AuxKernel() {}

  AuxKernel(const hephaestus::InputParameters &params){};

  virtual void
  Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
       hephaestus::DomainProperties &domain_properties){};

  virtual void Solve(double t = 0.0){};
};

class AuxKernels : public mfem::NamedFieldsMap<hephaestus::AuxKernel> {
public:
  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::DomainProperties &domain_properties);
  void Solve(double t = 0.0);
};

class VectorCoefficientAuxKernel : public AuxKernel {
public:
  VectorCoefficientAuxKernel(const hephaestus::InputParameters &params);

  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::DomainProperties &domain_properties) override;

  void Solve(double t = 0.0) override;

  std::string var_name;      // name of the variable
  std::string vec_coef_name; // name of the vector coefficient

  mfem::ParGridFunction *gf;
  mfem::VectorCoefficient *vec_coeff;
};

} // namespace hephaestus
