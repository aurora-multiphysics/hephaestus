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
  alphaCoef = new mfem::TransformedCoefficient(
      &oneCoef, domain_properties.scalar_property_map["magnetic_permeability"],
      fracFunc);
  betaCoef = domain_properties.scalar_property_map["electrical_conductivity"];
}

void ESolver::WriteConsoleSummary(double t, int it) {
  // Write a summary of the timestep to console.

  // Output Ohmic losses to console
  double el = this->ElectricLosses();
  if (myid_ == 0) {
    std::cout << std::fixed;
    std::cout << "step " << std::setw(6) << it << ",\tt = " << std::setw(6)
              << std::setprecision(3) << t
              << ",\tdot(E, J) = " << std::setprecision(8) << el << std::endl;
  }
}

double ESolver::ElectricLosses() const {
  double el = m1->InnerProduct(u_, u_);

  double global_el;
  MPI_Allreduce(&el, &global_el, 1, MPI_DOUBLE, MPI_SUM,
                m1->ParFESpace()->GetComm());

  return el;
}
} // namespace hephaestus
