#include "scaled_curl_vector_gridfunction_aux.hpp"

namespace hephaestus
{

ScaledCurlVectorGridFunctionAux::ScaledCurlVectorGridFunctionAux(
    const std::string & input_gf_name,
    const std::string & scaled_gf_name,
    const std::string & coef_name,
    const double & aConst,
    const hephaestus::InputParameters & solver_options)
  : ScaledVectorGridFunctionAux(input_gf_name, scaled_gf_name, coef_name, aConst, solver_options)
{
}

void
ScaledCurlVectorGridFunctionAux::BuildMixedBilinearForm()
{
  _a_mixed = std::make_unique<mfem::ParMixedBilinearForm>(_trial_fes, _test_fes);
  _a_mixed->AddDomainIntegrator(new mfem::MixedVectorCurlIntegrator(*_coef));
  _a_mixed->Assemble();
  _a_mixed->Finalize();
}

} // namespace hephaestus
