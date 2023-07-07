//* Solves:
//* ∇×(ν∇×A) + σdA/dt = Jᵉ
//*
//* in weak form
//* (ν∇×A, ∇×A') + (σdA/dt, A') - (Jᵉ, A') - <(ν∇×A)×n, A'>  = 0

//* where:
//* reluctivity ν = 1/μ
//* electrical_conductivity σ=1/ρ
//* Magnetic vector potential A
//* Electric field, E = ρJᵉ -dA/dt
//* Magnetic flux density, B = ∇×A
//* Magnetic field H = ν∇×A
//* Current density J = Jᵉ -σdA/dt

#include "a_formulation.hpp"

namespace hephaestus {

AFormulation::AFormulation() : HCurlFormulation() {
  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");
  h_curl_var_name = std::string("magnetic_vector_potential");
}

void AFormulation::RegisterAuxSolvers() {
  hephaestus::GridFunctions &variables = this->GetProblem()->gridfunctions;
  hephaestus::AuxSolvers &auxsolvers = this->GetProblem()->postprocessors;
  std::vector<std::string> aux_var_names;
  std::string b_field_name = "magnetic_flux_density";
  if (variables.Get(b_field_name) != NULL) {
    hephaestus::InputParameters b_field_aux_params;
    b_field_aux_params.SetParam("VariableName", h_curl_var_name);
    b_field_aux_params.SetParam("CurlVariableName", b_field_name);
    auxsolvers.Register("_magnetic_flux_density_aux",
                        new hephaestus::CurlAuxSolver(b_field_aux_params),
                        true);
  }
}

void AFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &domain_properties =
      this->GetProblem()->domain_properties;
  // std::cout << domain_properties.scalar_property_map.Has(
  //                  "electrical_conductivity")
  //           << std::endl;

  // mfem::IsoparametricTransformation T;
  // mfem::IntegrationPoint ip;

  // mfem::Coefficient *pw =
  //     domain_properties.scalar_property_map.Get("electrical_conductivity");
  // T.Attribute = 1;
  // std::cout << pw->Eval(T, ip) << std::endl;
  // T.Attribute = 2;
  // std::cout << pw->Eval(T, ip) << std::endl;

  // if (!domain_properties.scalar_property_map.Has("magnetic_permeability")) {
  //   domain_properties.scalar_property_map.Register(
  //       "magnetic_permeability",
  //       new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
  //           std::string("magnetic_permeability"))),
  //       true);
  // }
  // if (!domain_properties.scalar_property_map.Has("electrical_conductivity"))
  // {
  //   domain_properties.scalar_property_map.Register(
  //       "electrical_conductivity",
  //       new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
  //           std::string("electrical_conductivity"))),
  //       true);
  // }

  domain_properties.scalar_property_map.Register(
      alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map.Get("magnetic_permeability"),
          fracFunc),
      true);
}
} // namespace hephaestus
