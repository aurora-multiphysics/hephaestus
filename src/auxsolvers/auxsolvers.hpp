#pragma once
#include "auxsolver_base.hpp"
#include "coupled_coefficient_aux.hpp"
#include "curl_aux.hpp"
#include "helmholtz_projector.hpp"
#include "l2_error_vector_aux.hpp"
#include "scaled_curl_vector_gridfunction_aux.hpp"
#include "scaled_vector_gridfunction_aux.hpp"
#include "vector_coefficient_aux.hpp"
#include "vector_gridfunction_cross_product_aux.hpp"
#include "vector_gridfunction_dot_product_aux.hpp"

// Specify classes that perform auxiliary calculations on GridFunctions or
// Coefficients.
namespace hephaestus
{

class AuxSolvers : public hephaestus::NamedFieldsMap<hephaestus::AuxSolver>
{
private:
public:
  std::vector<hephaestus::AuxSolver *> _aux_queue;
  void Init(const hephaestus::GridFunctions & gridfunctions,
            hephaestus::Coefficients & coefficients);
  void Solve(double t = 0.0);
};

} // namespace hephaestus
