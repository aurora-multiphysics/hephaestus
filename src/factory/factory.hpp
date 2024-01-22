#pragma once
#include "../common/pfem_extras.hpp"
#include "a_formulation.hpp"
#include "av_formulation.hpp"
#include "complex_a_formulation.hpp"
#include "complex_e_formulation.hpp"
#include "e_formulation.hpp"
#include "eb_dual_formulation.hpp"
#include "h_formulation.hpp"
#include "inputs.hpp"
#include "magnetostatic_formulation.hpp"

namespace hephaestus
{

class Factory
{
public:
  static hephaestus::ProblemBuilder * CreateProblemBuilder(std::string & formulation_name);

  static hephaestus::TimeDomainEMFormulation *
  CreateTimeDomainEmFormulation(std::string & formulation);

  static hephaestus::FrequencyDomainEMFormulation *
  CreateFrequencyDomainEmFormulation(std::string & formulation);

  static mfem::ParFiniteElementSpace * CreateParFESpace(const hephaestus::InputParameters params,
                                                        mfem::ParMesh & pmesh);
};

} // namespace hephaestus
