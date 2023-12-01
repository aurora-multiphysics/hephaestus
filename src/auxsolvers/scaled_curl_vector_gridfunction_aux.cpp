#include "scaled_curl_vector_gridfunction_aux.hpp"

namespace hephaestus {

ScaledCurlVectorGridFunctionAux::ScaledCurlVectorGridFunctionAux(
    const std::string &input_gf_name, const std::string &scaled_gf_name,
    const std::string &coef_name, const double &aConst,
    const hephaestus::InputParameters &solver_options)
    : ScaledVectorGridFunctionAux(input_gf_name, scaled_gf_name, coef_name,
                                  aConst, solver_options) {}

void ScaledCurlVectorGridFunctionAux::buildMixedBilinearForm() {
  a_mixed = new mfem::ParMixedBilinearForm(trial_fes, test_fes);
  a_mixed->AddDomainIntegrator(new mfem::MixedVectorCurlIntegrator(*coef));
  a_mixed->Assemble();
  a_mixed->Finalize();
}

} // namespace hephaestus
