#include "auxkernels.hpp"

namespace hephaestus {

VectorCoefficientAuxKernel::VectorCoefficientAuxKernel(
    const std::string &var_name_, const std::string &vec_coef_name_)
    : var_name(var_name_), vec_coef_name(vec_coef_name_) {}

void VectorCoefficientAuxKernel::Init(
    const hephaestus::VariableMap &variables,
    hephaestus::DomainProperties &domain_properties) {
  gf = variables.Get(var_name);
  vec_coeff = domain_properties.vector_property_map[vec_coef_name];
}

void VectorCoefficientAuxKernel::Solve(double t) {
  vec_coeff->SetTime(t);
  gf->ProjectCoefficient(*vec_coeff);
}

} // namespace hephaestus
