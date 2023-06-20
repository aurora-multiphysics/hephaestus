#include "time_domain_problem_builder.hpp"

namespace hephaestus {

std::vector<mfem::ParGridFunction *>
TimeDomainProblemBuilder::RegisterTimeDerivatives(
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
