#include "e_solver.hpp"

namespace hephaestus {

ESolver::ESolver(mfem::ParMesh &pmesh, int order, hephaestus::BCMap &bc_map,
                 hephaestus::DomainProperties &domain_properties)
    : HCurlSolver(pmesh, order, bc_map, domain_properties) {}

void ESolver::SetVariableNames() {
  p_name = "electric_potential";
  p_display_name = "Scalar Potential (V)";

  u_name = "electric_field";
  u_display_name = "Electric Field (E)";

  v_name = "magnetic_flux_density";
  v_display_name = "Magnetic Flux Density (B)";
}

void ESolver::SetMaterialCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  sigmaCoef = domain_properties.scalar_property_map["electrical_conductivity"];
  muCoef = domain_properties.scalar_property_map["magnetic_permeability"];
  muInvCoef = new mfem::TransformedCoefficient(&oneCoef, muCoef, fracFunc);
  dtMuInvCoef = new mfem::TransformedCoefficient(&dtCoef, muInvCoef, prodFunc);
}

double ESolver::ElectricLosses() const {
  double el = m1->InnerProduct(e_, e_);

  double global_el;
  MPI_Allreduce(&el, &global_el, 1, MPI_DOUBLE, MPI_SUM,
                m1->ParFESpace()->GetComm());

  return el;
}
} // namespace hephaestus
