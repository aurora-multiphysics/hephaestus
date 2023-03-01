#include "factory.hpp"

namespace hephaestus {

hephaestus::TransientFormulation *
Factory::createTransientFormulation(std::string &formulation) {
  if (formulation == "EBForm") {
    return new hephaestus::EBDualFormulation();
  } else if (formulation == "HJForm") {
    return new hephaestus::HJDualFormulation();
  } else if (formulation == "HForm") {
    return new hephaestus::HFormulation();
  } else if (formulation == "AForm") {
    return new hephaestus::AFormulation();
  } else if (formulation == "EForm") {
    return new hephaestus::EFormulation();
  } else if (formulation == "AVForm") {
    return new hephaestus::AVFormulation();
  } else {
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
