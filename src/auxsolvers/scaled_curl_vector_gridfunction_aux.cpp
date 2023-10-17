#include "scaled_curl_vector_gridfunction_aux.hpp"

namespace hephaestus {

ScaledCurlVectorGridFunctionAux::ScaledCurlVectorGridFunctionAux(
    const hephaestus::InputParameters &params)
    : ScaledVectorGridFunctionAux(params) {}

void ScaledCurlVectorGridFunctionAux::buildMixedBilinearForm() {
  a_mixed = new mfem::ParMixedBilinearForm(trial_fes, test_fes);
  a_mixed->AddDomainIntegrator(new mfem::MixedVectorCurlIntegrator(*coef));
  a_mixed->Assemble();
}

} // namespace hephaestus
