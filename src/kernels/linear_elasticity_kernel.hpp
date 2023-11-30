#pragma once
#include "kernel_base.hpp"

namespace hephaestus {

/*
Implements the provided mfem linear elasticity integrator into a hephaestus kernel.
a(u,v) = (λ div(u), div(v)) + (2 μ e(u), e(v)),
where e(v) = (1/2) (grad(v) + grad(v)^T).
*/
class LinearElasticityKernel : public Kernel<mfem::ParBilinearForm> {
public:
  LinearElasticityKernel(const hephaestus::InputParameters &params);
  virtual void Init(hephaestus::GridFunctions &gridfunctions,
                    const hephaestus::FESpaces &fespaces,
                    hephaestus::BCMap &bc_map,
                    hephaestus::Coefficients &coefficients) override;
  virtual void Apply(mfem::ParBilinearForm *blf) override;
  std::string lame_parameter_name, shear_modulus_name;
  mfem::Coefficient *lame_parameter, *shear_modulus;

};

}; // namespace hephaestus
