#pragma once
#include "../common/pfem_extras.hpp"
#include "aform_solver.hpp"
#include "av_solver.hpp"
#include "complex_a_form.hpp"
#include "complex_e_form.hpp"
#include "eb_dual_solver.hpp"
#include "eform_solver.hpp"
#include "hform_solver.hpp"
#include "hj_dual_solver.hpp"
#include "inputs.hpp"

namespace hephaestus {

class Factory {
public:
  static hephaestus::TimeDomainFormulation *
  createTimeDomainFormulation(std::string &formulation);

  static hephaestus::FrequencyDomainFormulation *
  createFrequencyDomainFormulation(std::string &formulation);

  static mfem::ParFiniteElementSpace *
  createParFESpace(const hephaestus::InputParameters params,
                   mfem::ParMesh &pmesh);
};

} // namespace hephaestus
