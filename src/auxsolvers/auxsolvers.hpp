#pragma once
#include "coefficient_aux.hpp"
#include "coupled_coefficient_aux.hpp"
#include "curl_aux.hpp"
#include "joule_heating_aux.hpp"
#include "l2_error_vector_aux.hpp"
#include "scaled_curl_vector_gridfunction_aux.hpp"
#include "scaled_gridfunction_aux.hpp"
#include "scaled_vector_gridfunction_aux.hpp"
#include "vector_coefficient_aux.hpp"

// Specify classes that perform auxiliary calculations on GridFunctions or
// Coefficients.
namespace hephaestus {

class AuxSolvers : public mfem::NamedFieldsMap<hephaestus::AuxSolver> {
private:
public:
  std::vector<hephaestus::AuxSolver *> aux_queue;
  void Init(const hephaestus::GridFunctions &gridfunctions,
            hephaestus::Coefficients &coefficients);
  void Solve(double t = 0.0);
};

} // namespace hephaestus
