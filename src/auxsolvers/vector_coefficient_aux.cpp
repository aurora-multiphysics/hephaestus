#include "vector_coefficient_aux.hpp"

namespace hephaestus {

VectorCoefficientAuxSolver::VectorCoefficientAuxSolver(
    const hephaestus::InputParameters &params)
    : AuxSolver(), var_name(params.GetParam<std::string>("VariableName")),
      vec_coef_name(params.GetParam<std::string>("VectorCoefficientName")) {}

void VectorCoefficientAuxSolver::Init(
    const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::DomainProperties &domain_properties) {
  gf = variables.Get(var_name);
  vec_coeff = domain_properties.vector_property_map.Get(vec_coef_name);
}

void VectorCoefficientAuxSolver::Solve(double t) {
  vec_coeff->SetTime(t);
  gf->ProjectCoefficient(*vec_coeff);
}

} // namespace hephaestus
