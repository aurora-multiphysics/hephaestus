#include "factory.hpp"

namespace hephaestus {

hephaestus::TransientFormulation *Factory::createTransientFormulation(
    std::string &formulation, mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options) {
  if (formulation == "EBForm") {
    return new hephaestus::EBDualSolver(pmesh, order, fespaces, variables,
                                        bc_map, domain_properties, sources,
                                        solver_options);
  } else if (formulation == "HJForm") {
    return new hephaestus::HJDualSolver(pmesh, order, fespaces, variables,
                                        bc_map, domain_properties, sources,
                                        solver_options);
  } else if (formulation == "HForm") {
    return new hephaestus::HFormSolver(pmesh, order, fespaces, variables,
                                       bc_map, domain_properties, sources,
                                       solver_options);
  } else if (formulation == "AForm") {
    return new hephaestus::AFormSolver(pmesh, order, fespaces, variables,
                                       bc_map, domain_properties, sources,
                                       solver_options);
  } else if (formulation == "EForm") {
    return new hephaestus::EFormSolver(pmesh, order, fespaces, variables,
                                       bc_map, domain_properties, sources,
                                       solver_options);
  } else if (formulation == "AVForm") {
    return new hephaestus::AVSolver(pmesh, order, fespaces, variables, bc_map,
                                    domain_properties, sources, solver_options);
  } else if (formulation == "LinearElasticSolver") {
    return new hephaestus::LinearElasticSolver(pmesh, order, fespaces, variables, bc_map,
                                    domain_properties, sources, solver_options);
  }
   else {
    std::cout << "Formulation name " << formulation << " not recognised. \n";
  }
  return nullptr;
}

mfem::ParFiniteElementSpace *
Factory::createParFESpace(hephaestus::InputParameters params,
                          mfem::ParMesh &pmesh) {
  std::string FEType(params.GetParam<std::string>("FESpaceType"));
  int order(params.GetParam<int>("order"));
  int components(params.GetParam<int>(
      "components")); // spatial dimension of mesh. Use
                      // FiniteElementCollection::New instead
  if (FEType == "H1") {
    return new mfem::common::H1_ParFESpace(&pmesh, order, components);
  } else if (FEType == "ND") {
    return new mfem::common::ND_ParFESpace(&pmesh, order, components);
  } else if (FEType == "RT") {
    return new mfem::common::RT_ParFESpace(&pmesh, order, components);
  } else if (FEType == "L2") {
    return new mfem::common::L2_ParFESpace(&pmesh, order, components);
  } else {
    std::cout << "FESpaceType " << FEType << " not recognised. \n";
  }
  return nullptr;
}

} // namespace hephaestus
