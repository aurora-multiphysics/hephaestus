#pragma once
#include "coefficient_aux.hpp"
#include "coupled_coefficient_aux.hpp"
#include "curl_aux.hpp"
#include "l2_error_vector_aux.hpp"
#include "vector_coefficient_aux.hpp"

// Specify classes that perform auxiliary calculations on GridFunctions or
// Coefficients.
namespace hephaestus {

class AuxSolvers : public mfem::NamedFieldsMap<hephaestus::AuxSolver> {
private:
public:
  std::vector<hephaestus::AuxSolver *> aux_queue;
  void Init(const mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
            hephaestus::DomainProperties &domain_properties);
  void Solve(double t = 0.0);
};

} // namespace hephaestus
