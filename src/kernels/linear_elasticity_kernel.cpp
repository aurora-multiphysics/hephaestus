#include "linear_elasticity_kernel.hpp"

namespace hephaestus {

LinearElasticityKernel::LinearElasticityKernel(
    const hephaestus::InputParameters &params)
    : Kernel(params),
      lame_paramter_name(params.GetParam<std::string>("LameParameterCoefName")),
      shear_modulus_name(params.GetParam<std::string>("ShearModulusCoefName"))
    {}

void LinearElasticityKernel::Init(hephaestus::GridFunctions &gridfunctions,
                                     const hephaestus::FESpaces &fespaces,
                                     hephaestus::BCMap &bc_map,
                                     hephaestus::Coefficients &coefficients) {

  lame_parameter = coefficients.scalars.Get(lame_paramter_name);
  shear_modulus = coefficients.scalars.Get(shear_modulus_name);
}

void LinearElasticityKernel::Apply(mfem::ParBilinearForm *blf) {
  blf->AddDomainIntegrator(new mfem::ElasticityIntegrator(*lame_parameter, *shear_modulus));
}

} // namespace hephaestus
