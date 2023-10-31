#include "open_coil.hpp"

namespace hephaestus {

    OpenCoilSolver::OpenCoilSolver(
    const hephaestus::InputParameters &params,
    const std::vector<hephaestus::Subdomain> &coil_dom,
    const std::pair<int,int> electrodes, const int order)
    : hcurl_fespace_name_(params.GetParam<std::string>("HCurlFESpaceName")),
      h1_fespace_name_(params.GetParam<std::string>("H1FESpaceName")),
      J_gf_name_(params.GetParam<std::string>("JGridFunctionName")),
      V_gf_name_(params.GetParam<std::string>("VGridFunctionName")),
      I_coef_name_(params.GetParam<std::string>("IFuncCoefName")),
      coil_domains_(coil_dom), order_(order), elec_attrs_(electrodes),coef1_(nullptr), coef0_(nullptr),
      mesh_parent_(nullptr), J_parent_(nullptr), V_parent_(nullptr), HCurlFESpace_parent_(nullptr),
      H1FESpace_parent_(nullptr) {

  coef1_ = new mfem::ConstantCoefficient(1.0);
  coef0_ = new mfem::ConstantCoefficient(0.0);
}


OpenCoilSolver::~OpenCoilSolver() {
}


void OpenCoilSolver::Init(hephaestus::GridFunctions &gridfunctions,
                            const hephaestus::FESpaces &fespaces,
                            hephaestus::BCMap &bc_map,
                            hephaestus::Coefficients &coefficients) {
}

void OpenCoilSolver::Apply(mfem::ParLinearForm *lf) {

}

void OpenCoilSolver::SubtractSource(mfem::ParGridFunction *gf) {}


} // namespace hephaestus