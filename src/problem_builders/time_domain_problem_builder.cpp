#include "time_domain_problem_builder.hpp"

namespace hephaestus {
TransientProblem::TransientProblem(const hephaestus::InputParameters &params)
    : Problem(params) {}

std::unique_ptr<hephaestus::TimeDependentEquationSystem>
TransientProblemBuilder::CreateTimeDependentEquationSystem() const {
  hephaestus::InputParameters params;
  return std::make_unique<hephaestus::TimeDependentEquationSystem>(params);
};

std::unique_ptr<hephaestus::TimeDomainEquationSystemOperator>
TransientProblemBuilder::CreateTimeDomainEquationSystemOperator(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources,
    hephaestus::InputParameters &solver_options) const {
  return std::make_unique<hephaestus::TimeDomainEquationSystemOperator>(
      pmesh, fespaces, variables, bc_map, domain_properties, sources,
      solver_options);
};

std::vector<mfem::ParGridFunction *>
TransientProblemBuilder::RegisterTimeDerivatives(
    std::vector<std::string> gridfunction_names,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &gridfunctions) {
  std::vector<mfem::ParGridFunction *> time_derivatives;

  for (auto &gridfunction_name : gridfunction_names) {
    gridfunctions.Register(
        GetTimeDerivativeName(gridfunction_name),
        new mfem::ParGridFunction(
            gridfunctions.Get(gridfunction_name)->ParFESpace()),
        true);
    time_derivatives.push_back(
        gridfunctions.Get(GetTimeDerivativeName(gridfunction_name)));
  }
  return time_derivatives;
}

} // namespace hephaestus
