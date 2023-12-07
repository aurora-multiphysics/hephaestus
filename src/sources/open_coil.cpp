#include "open_coil.hpp"
#include "utils.hpp"

namespace hephaestus {

double highV(const mfem::Vector &x, double t) { return 1.0; }
double lowV(const mfem::Vector &x, double t) { return 0.0; }

OpenCoilSolver::OpenCoilSolver(const hephaestus::InputParameters &params,
                               const mfem::Array<int> &coil_dom,
                               const std::pair<int, int> electrodes)
    : J_gf_name_(params.GetParam<std::string>("SourceName")),
      V_gf_name_(params.GetParam<std::string>("PotentialName")),
      I_coef_name_(params.GetParam<std::string>("IFuncCoefName")),
      coil_domains_(coil_dom), elec_attrs_(electrodes), coef1_(1.0),
      mesh_parent_(nullptr), mesh_(nullptr), H1FESpace_(nullptr),
      HCurlFESpace_(nullptr), J_parent_(nullptr), Jt_parent_(nullptr),
      V_parent_(nullptr), Vt_parent_(nullptr), J_(nullptr), V_(nullptr),
      m1_(nullptr), final_lf_(nullptr), high_src_(highV), low_src_(lowV),
      high_terminal_(1), low_terminal_(1) {

  hephaestus::InputParameters default_pars;
  default_pars.SetParam("Tolerance", float(1.0e-20));
  default_pars.SetParam("AbsTolerance", float(1.0e-20));
  default_pars.SetParam("MaxIter", (unsigned int)1000);
  default_pars.SetParam("PrintLevel", 1);

  solver_options_ = params.GetOptionalParam<hephaestus::InputParameters>(
      "SolverOptions", default_pars);

  ref_face_ = elec_attrs_.first;
}

OpenCoilSolver::~OpenCoilSolver() {

  ifDelete(mesh_);
  ifDelete(m1_);
  ifDelete(H1FESpace_);
  ifDelete(HCurlFESpace_);
  ifDelete(J_);
  ifDelete(V_);
  ifDelete(Jt_parent_);
  ifDelete(Vt_parent_);
  ifDelete(final_lf_);
}

void OpenCoilSolver::Init(hephaestus::GridFunctions &gridfunctions,
                          const hephaestus::FESpaces &fespaces,
                          hephaestus::BCMap &bc_map,
                          hephaestus::Coefficients &coefficients) {

  Itotal_ = coefficients.scalars.Get(I_coef_name_);
  if (Itotal_ == nullptr) {
    std::cout << I_coef_name_ + " not found in coefficients when "
                                "creating OpenCoilSolver. "
                                "Assuming unit current.\n";
    Itotal_ = new mfem::ConstantCoefficient(1.0);
  }

  J_parent_ = gridfunctions.Get(J_gf_name_);
  if (J_parent_ == nullptr) {
    const std::string error_message = J_gf_name_ +
                                      " not found in gridfunctions when "
                                      "creating OpenCoilSolver\n";
    mfem::mfem_error(error_message.c_str());
  } else if (J_parent_->ParFESpace()->FEColl()->GetContType() !=
             mfem::FiniteElementCollection::TANGENTIAL) {
    mfem::mfem_error("J GridFunction must be of HCurl type.");
  }
  order_hcurl_ = J_parent_->ParFESpace()->FEColl()->GetOrder();

  V_parent_ = gridfunctions.Get(V_gf_name_);
  if (V_parent_ == nullptr) {
    std::cout << V_gf_name_ + " not found in gridfunctions when "
                              "creating OpenCoilSolver.\n";
    order_h1_ = order_hcurl_;
  } else if (V_parent_->ParFESpace()->FEColl()->GetContType() !=
             mfem::FiniteElementCollection::CONTINUOUS) {
    mfem::mfem_error("V GridFunction must be of H1 type.");
  } else {
    order_h1_ = V_parent_->ParFESpace()->FEColl()->GetOrder();
    Vt_parent_ = new mfem::ParGridFunction(*V_parent_);
  }

  mesh_parent_ = J_parent_->ParFESpace()->GetParMesh();

  initChildMesh();
  makeFESpaces();
  makeGridFunctions();
  setBCs();
  SPSCurrent();
}

void OpenCoilSolver::Apply(mfem::ParLinearForm *lf) {

  // The transformation and integration points themselves are not relevant, it's
  // just so we can call Eval
  mfem::ElementTransformation *Tr = mesh_parent_->GetElementTransformation(0);
  const mfem::IntegrationPoint &ip =
      mfem::IntRules.Get(J_parent_->ParFESpace()->GetFE(0)->GetGeomType(), 1)
          .IntPoint(0);

  double I = Itotal_->Eval(*Tr, ip);

  *J_parent_ = 0.0;
  J_parent_->Add(I, *Jt_parent_);
  if (V_parent_ != nullptr) {
    *V_parent_ = 0.0;
    V_parent_->Add(I, *Vt_parent_);
  }

  lf->Add(I, *final_lf_);
}

void OpenCoilSolver::SubtractSource(mfem::ParGridFunction *gf) {}

void OpenCoilSolver::initChildMesh() {

  if (mesh_ == nullptr)
    mesh_ = new mfem::ParSubMesh(
        mfem::ParSubMesh::CreateFromDomain(*mesh_parent_, coil_domains_));
}

void OpenCoilSolver::makeFESpaces() {

  if (H1FESpace_ == nullptr)
    H1FESpace_ = new mfem::ParFiniteElementSpace(
        mesh_, new mfem::H1_FECollection(order_h1_, mesh_->Dimension()));

  if (HCurlFESpace_ == nullptr)
    HCurlFESpace_ = new mfem::ParFiniteElementSpace(
        mesh_, new mfem::ND_FECollection(order_hcurl_, mesh_->Dimension()));
}

void OpenCoilSolver::makeGridFunctions() {

  if (V_ == nullptr)
    V_ = new mfem::ParGridFunction(H1FESpace_);

  if (J_ == nullptr)
    J_ = new mfem::ParGridFunction(HCurlFESpace_);
  
  if (Jt_parent_ == nullptr)
    Jt_parent_ = new mfem::ParGridFunction(*J_parent_);

  *V_ = 0.0;
  *J_ = 0.0;
  *Jt_parent_ = 0.0;
}

void OpenCoilSolver::setBCs() {

  high_terminal_[0] = elec_attrs_.first;
  low_terminal_[0] = elec_attrs_.second;
}

void OpenCoilSolver::SPSCurrent() {

  bc_maps.Register("high_potential",
                   new hephaestus::ScalarDirichletBC(
                       std::string("V"), high_terminal_, &high_src_),
                   true);

  bc_maps.Register("low_potential",
                   new hephaestus::ScalarDirichletBC(std::string("V"),
                                                     low_terminal_, &low_src_),
                   true);

  hephaestus::FESpaces fespaces;
  fespaces.Register(std::string("HCurl"), HCurlFESpace_, true);
  fespaces.Register(std::string("H1"), H1FESpace_, true);

  hephaestus::GridFunctions gridfunctions;
  gridfunctions.Register(std::string("source"), J_, true);
  gridfunctions.Register(std::string("V"), V_, true);

  hephaestus::InputParameters sps_params;
  sps_params.SetParam("SourceName", std::string("source"));
  sps_params.SetParam("PotentialName", std::string("V"));
  sps_params.SetParam("HCurlFESpaceName", std::string("HCurl"));
  sps_params.SetParam("H1FESpaceName", std::string("H1"));
  sps_params.SetParam("SolverOptions", solver_options_);
  sps_params.SetParam("ConductivityCoefName",
                      std::string("magnetic_permeability"));

  hephaestus::Coefficients coefs;
  coefs.scalars.Register("magnetic_permeability", &coef1_, false);

  hephaestus::ScalarPotentialSource sps(sps_params);
  sps.Init(gridfunctions, fespaces, bc_maps, coefs);

  mfem::ParLinearForm dummy(HCurlFESpace_);
  sps.Apply(&dummy);

  // Normalise the current through the wedges and use them as a reference
  double flux = calcFlux(J_, ref_face_);
  *J_ /= abs(flux);
  if (V_)
    *V_ /= abs(flux);

  mesh_->Transfer(*J_, *Jt_parent_);
  if (V_parent_)
    mesh_->Transfer(*V_,*Vt_parent_);

  buildM1();

  final_lf_ = new mfem::ParLinearForm(Jt_parent_->ParFESpace());
  *final_lf_ = 0.0;
  m1_->AddMult(*Jt_parent_, *final_lf_, 1.0);
}

void OpenCoilSolver::buildM1() {

  if (m1_ == nullptr) {

    m1_ = new mfem::ParBilinearForm(J_parent_->ParFESpace());
    hephaestus::attrToMarker(coil_domains_, coil_markers_,
                             mesh_parent_->attributes.Max());
    m1_->AddDomainIntegrator(
        new mfem::VectorFEMassIntegrator(new mfem::ConstantCoefficient(1.0)),
        coil_markers_);
    m1_->Assemble();
    m1_->Finalize();
  }
}

void OpenCoilSolver::setRefFace(const int face) { ref_face_ = face; }

} // namespace hephaestus