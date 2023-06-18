#pragma once
#include "factory.hpp"
#include "frequency_domain_problem_builder.hpp"
#include "time_domain_problem_builder.hpp"

namespace hephaestus {

static hephaestus::ProblemBuilder *
createProblemBuilder(std::string &formulation_name) {
  if (formulation_name == "EBForm" || formulation_name == "HJForm" ||
      formulation_name == "HForm" || formulation_name == "AForm" ||
      formulation_name == "EForm" || formulation_name == "AVForm" ||
      formulation_name == "Custom") {
    hephaestus::TransientProblemBuilder *builder =
        new hephaestus::TransientProblemBuilder();
    // hephaestus::TransientFormulation *_formulation =
    //     hephaestus::Factory::createTransientFormulation(formulation_name);
    // builder->SetFormulation(_formulation);
    return builder;
  } else if (formulation_name == "ComplexEForm" ||
             formulation_name == "ComplexAForm") {
    hephaestus::FrequencyDomainProblemBuilder *builder =
        new hephaestus::FrequencyDomainProblemBuilder();
    ;
    // hephaestus::SteadyFormulation *_formulation =
    //     hephaestus::Factory::createSteadyFormulation();
    // builder->SetFormulation(_formulation);
    return builder;

    // return new hephaestus::FrequencyDomainProblemBuilder();
  } else {
    std::cout << "Formulation name " << formulation_name
              << " not recognised. \n";
  }
  return nullptr;
};

} // namespace hephaestus
