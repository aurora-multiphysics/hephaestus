#include "factory.hpp"

namespace hephaestus {

hephaestus::ProblemBuilder *
Factory::createProblemBuilder(std::string &formulation_name) {
  if (formulation_name == "EBForm") {
    return new hephaestus::EBDualFormulation(
        "magnetic_reluctivity", "magnetic_permeability",
        "electrical_conductivity", "electric_field", "magnetic_flux_density");
  } else if (formulation_name == "HJForm") {
    return new hephaestus::HJDualFormulation(
        "electrical_resistivity", "electrical_conductivity",
        "magnetic_permeability", "magnetic_field", "current_density");
  } else if (formulation_name == "HForm") {
    return new hephaestus::HFormulation(
        "electrical_resistivity", "electrical_conductivity",
        "magnetic_permeability", "magnetic_field");
  } else if (formulation_name == "AForm") {
    return new hephaestus::AFormulation(
        "magnetic_reluctivity", "magnetic_permeability",
        "electrical_conductivity", "magnetic_vector_potential");
  } else if (formulation_name == "EForm") {
    return new hephaestus::EFormulation(
        "magnetic_reluctivity", "magnetic_permeability",
        "electrical_conductivity", "electric_field");
  } else if (formulation_name == "AVForm") {
    return new hephaestus::AVFormulation(
        "magnetic_reluctivity", "magnetic_permeability",
        "electrical_conductivity", "magnetic_vector_potential",
        "electric_potential");
  } else if (formulation_name == "ComplexEForm") {
    return new hephaestus::ComplexEFormulation(
        "magnetic_reluctivity", "electrical_conductivity",
        "dielectric_permittivity", "frequency", "electric_field");
  } else if (formulation_name == "ComplexAForm") {
    return new hephaestus::ComplexAFormulation(
        "magnetic_reluctivity", "electrical_conductivity",
        "dielectric_permittivity", "frequency", "magnetic_vector_potential");
  } else if (formulation_name == "Custom") {
    return new hephaestus::TimeDomainEMFormulation();
  } else {
    MFEM_WARNING("Formulation name " << formulation_name << " not recognised.");
    return nullptr;
  }
};

hephaestus::FrequencyDomainEMFormulation *
Factory::createFrequencyDomainEMFormulation(std::string &formulation) {
  if (formulation == "ComplexEForm") {
    return new hephaestus::ComplexEFormulation(
        "magnetic_reluctivity", "electrical_conductivity",
        "dielectric_permittivity", "frequency", "electric_field");
  } else if (formulation == "ComplexAForm") {
    return new hephaestus::ComplexAFormulation(
        "magnetic_reluctivity", "electrical_conductivity",
        "dielectric_permittivity", "frequency", "magnetic_vector_potential");
  } else {
    MFEM_WARNING("Steady formulation name " << formulation
                                            << " not recognised.");
  }
  return nullptr;
}

hephaestus::TimeDomainEMFormulation *
Factory::createTimeDomainEMFormulation(std::string &formulation) {
  if (formulation == "EBForm") {
    return new hephaestus::EBDualFormulation(
        "magnetic_reluctivity", "magnetic_permeability",
        "electrical_conductivity", "electric_field", "magnetic_flux_density");
  } else if (formulation == "HJForm") {
    return new hephaestus::HJDualFormulation(
        "electrical_resistivity", "electrical_conductivity",
        "magnetic_permeability", "magnetic_field", "current_density");
  } else if (formulation == "HForm") {
    return new hephaestus::HFormulation(
        "electrical_resistivity", "electrical_conductivity",
        "magnetic_permeability", "magnetic_field");
  } else if (formulation == "AForm") {
    return new hephaestus::AFormulation(
        "magnetic_reluctivity", "magnetic_permeability",
        "electrical_conductivity", "magnetic_vector_potential");
  } else if (formulation == "EForm") {
    return new hephaestus::EFormulation(
        "magnetic_reluctivity", "magnetic_permeability",
        "electrical_conductivity", "electric_field");
  } else if (formulation == "AVForm") {
    return new hephaestus::AVFormulation(
        "magnetic_reluctivity", "magnetic_permeability",
        "electrical_conductivity", "magnetic_vector_potential",
        "electric_potential");
  } else if (formulation == "Custom") {
    return new hephaestus::TimeDomainEMFormulation();
  } else {
    MFEM_WARNING("Formulation name " << formulation << " not recognised.");
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
    MFEM_WARNING("FESpaceType " << FEType << " not recognised.");
  }
  return nullptr;
}

} // namespace hephaestus
