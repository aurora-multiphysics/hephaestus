#include "factory.hpp"

namespace hephaestus {

hephaestus::ProblemBuilder *
Factory::createProblemBuilder(std::string &formulation_name) {
  if (formulation_name == "EBForm") {
    return new hephaestus::EBDualFormulation();
  } else if (formulation_name == "HJForm") {
    return new hephaestus::HJDualFormulation();
  } else if (formulation_name == "HForm") {
    return new hephaestus::HFormulation();
  } else if (formulation_name == "AForm") {
    return new hephaestus::AFormulation();
  } else if (formulation_name == "EForm") {
    return new hephaestus::EFormulation();
  } else if (formulation_name == "AVForm") {
    return new hephaestus::AVFormulation();
  } else if (formulation_name == "ComplexEForm") {
    return new hephaestus::ComplexEFormulation();
  } else if (formulation_name == "ComplexAForm") {
    return new hephaestus::ComplexAFormulation();
  } else if (formulation_name == "Custom") {
    return new hephaestus::TimeDomainFormulation();
  } else {
    MFEM_WARNING("Formulation name " << formulation_name << " not recognised.");
    return nullptr;
  }
};

hephaestus::FrequencyDomainFormulation *
Factory::createFrequencyDomainFormulation(std::string &formulation) {
  if (formulation == "ComplexEForm") {
    return new hephaestus::ComplexEFormulation();
  } else if (formulation == "ComplexAForm") {
    return new hephaestus::ComplexAFormulation();
  } else {
    std::cout << "Steady formulation name " << formulation
              << " not recognised. \n";
  }
  return nullptr;
}

hephaestus::TimeDomainFormulation *
Factory::createTimeDomainFormulation(std::string &formulation) {
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
  } else if (formulation == "Custom") {
    return new hephaestus::TimeDomainFormulation();
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
