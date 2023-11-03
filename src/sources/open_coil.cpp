#include "open_coil.hpp"

namespace hephaestus {

double highV(const mfem::Vector &x, double t) { return 1.0; }
double lowV(const mfem::Vector &x, double t) { return 0.0; }

///// THESE FUNCTIONS WILL EVENTUALLY GO INTO A UTILS FILE ///////////

double calcFlux(mfem::GridFunction *v_field, int face_attr) {

  double flux = 0.0;
  double area = 0.0;

  mfem::FiniteElementSpace *FES = v_field->FESpace();
  mfem::Mesh *mesh = FES->GetMesh();

  mfem::Vector local_dofs, normal_vec;
  mfem::DenseMatrix dshape;
  mfem::Array<int> dof_ids;

  for (int i = 0; i < mesh->GetNBE(); i++) {

    if (mesh->GetBdrAttribute(i) != face_attr)
      continue;

    mfem::FaceElementTransformations *FTr =
        mesh->GetFaceElementTransformations(mesh->GetBdrFace(i));
    if (FTr == nullptr)
      continue;

    const mfem::FiniteElement &elem = *FES->GetFE(FTr->Elem1No);
    const int int_order = 2 * elem.GetOrder() + 3;
    const mfem::IntegrationRule &ir =
        mfem::IntRules.Get(FTr->FaceGeom, int_order);

    FES->GetElementDofs(FTr->Elem1No, dof_ids);
    v_field->GetSubVector(dof_ids, local_dofs);
    const int space_dim = FTr->Face->GetSpaceDim();
    normal_vec.SetSize(space_dim);
    dshape.SetSize(elem.GetDof(), space_dim);

    for (int j = 0; j < ir.GetNPoints(); j++) {

      const mfem::IntegrationPoint &ip = ir.IntPoint(j);
      mfem::IntegrationPoint eip;
      FTr->Loc1.Transform(ip, eip);
      FTr->Face->SetIntPoint(&ip);
      double face_weight = FTr->Face->Weight();
      double val = 0.0;
      FTr->Elem1->SetIntPoint(&eip);
      elem.CalcVShape(*FTr->Elem1, dshape);
      mfem::CalcOrtho(FTr->Face->Jacobian(), normal_vec);
      val += dshape.InnerProduct(normal_vec, local_dofs) / face_weight;

      // Measure the area of the boundary
      area += ip.weight * face_weight;

      // Integrate alpha * n.Grad(x) + beta * x
      flux += val * ip.weight * face_weight;
    }
  }

  double total_flux;
  MPI_Allreduce(&flux, &total_flux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return total_flux;
}

void SubdomainToArray(const std::vector<hephaestus::Subdomain> &sd,
                      mfem::Array<int> &arr) {
  arr.DeleteAll();
  for (auto s : sd)
    arr.Append(s.id);
}

void SubdomainToArray(const hephaestus::Subdomain &sd, mfem::Array<int> &arr) {
  arr.DeleteAll();
  arr.Append(sd.id);
}

template <typename T> void ifDelete(T *ptr) {
  if (ptr != nullptr)
    delete ptr;
}

void inheritBdrAttributes(const mfem::ParMesh *parent_mesh,
                          mfem::ParSubMesh *child_mesh) {

  int face, ori, att;
  auto map = child_mesh->GetParentToSubMeshFaceIDMap();

  for (int bdr = 0; bdr < parent_mesh->GetNBE(); ++bdr) {

    parent_mesh->GetBdrElementFace(bdr, &face, &ori);
    if (map[face] != -1) {
      att = parent_mesh->GetBdrAttribute(bdr);
      auto *new_elem = child_mesh->GetFace(map[face])->Duplicate(child_mesh);
      new_elem->SetAttribute(att);
      child_mesh->AddBdrElement(new_elem);
    }
  }

  child_mesh->FinalizeTopology();
  child_mesh->Finalize();
  child_mesh->SetAttributes();
}

void cleanDivergence(hephaestus::GridFunctions *gridfunctions,
                     std::string J_name, std::string V_name,
                     hephaestus::BCMap *bc_map) {

  hephaestus::InputParameters pars;
  hephaestus::FESpaces fes;

  pars.SetParam("VectorGridFunctionName", J_name);
  pars.SetParam("ScalarGridFunctionName", V_name);
  hephaestus::HelmholtzProjector projector(pars);
  projector.Project(*gridfunctions, fes, *bc_map);
}

/////////////////////////////////////////////////////////////////////

OpenCoilSolver::OpenCoilSolver(const hephaestus::InputParameters &params,
                               const mfem::Array<int> &coil_dom,
                               const std::pair<int, int> electrodes)
    : J_gf_name_(params.GetParam<std::string>("SourceName")),
      V_gf_name_(params.GetParam<std::string>("PotentialName")),
      I_coef_name_(params.GetParam<std::string>("IFuncCoefName")),
      coil_domains_(coil_dom), elec_attrs_(electrodes), coef1_(1.0),
      mesh_parent_(nullptr), J_parent_(nullptr), V_parent_(nullptr),
      J_(nullptr), V_(nullptr), high_src_(highV), low_src_(lowV) {

  ref_face_ = elec_attrs_.first;
}

OpenCoilSolver::~OpenCoilSolver() {

  ifDelete(mesh_);
  ifDelete(H1_Collection_);
  ifDelete(HCurl_Collection_);
  ifDelete(H1FESpace_);
  ifDelete(HCurlFESpace_);

  ifDelete(J_);
  ifDelete(V_);

  ifDelete(bc_maps_);
  ifDelete(high_DBC_);
  ifDelete(low_DBC_);
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
  *J_ *= I;
  mesh_->Transfer(*J_, *J_parent_);
  *J_ /= I;

  if (V_parent_ != nullptr) {
    *V_ *= I;
    mesh_->Transfer(*V_, *V_parent_);
    *V_ /= I;
  }

  lf->Add(1.0, *J_parent_);
}

void OpenCoilSolver::SubtractSource(mfem::ParGridFunction *gf) {}

void OpenCoilSolver::initChildMesh() {

  mesh_ = new mfem::ParSubMesh(
      mfem::ParSubMesh::CreateFromDomain(*mesh_parent_, coil_domains_));

  inheritBdrAttributes(mesh_parent_, mesh_);
}

void OpenCoilSolver::makeFESpaces() {

  H1_Collection_ = new mfem::H1_FECollection(order_h1_, mesh_->Dimension());
  HCurl_Collection_ =
      new mfem::ND_FECollection(order_hcurl_, mesh_->Dimension());
  H1FESpace_ = new mfem::ParFiniteElementSpace(mesh_, H1_Collection_);
  HCurlFESpace_ = new mfem::ParFiniteElementSpace(mesh_, HCurl_Collection_);
}

void OpenCoilSolver::makeGridFunctions() {

  if (V_ == nullptr)
    V_ = new mfem::ParGridFunction(H1FESpace_);

  if (J_ == nullptr)
    J_ = new mfem::ParGridFunction(HCurlFESpace_);

  *V_ = 0.0;
  *J_ = 0.0;
}

void OpenCoilSolver::setBCs() {

  if (high_terminal_.Size() == 0)
    high_terminal_.Append(elec_attrs_.first);
  if (low_terminal_.Size() == 0)
    low_terminal_.Append(elec_attrs_.second);

  high_DBC_ = new hephaestus::ScalarDirichletBC(std::string("V"),
                                                  high_terminal_, &high_src_);
  low_DBC_ = new hephaestus::ScalarDirichletBC(std::string("V"),
                                                 low_terminal_, &low_src_);

  bc_maps_ = new hephaestus::BCMap;
  bc_maps_->Register("high_potential", high_DBC_, true);
  bc_maps_->Register("low_potential", low_DBC_, true);
}

void OpenCoilSolver::SPSCurrent() {

  fespaces_ = new hephaestus::FESpaces;
  fespaces_->Register(std::string("HCurl"), HCurlFESpace_, true);
  fespaces_->Register(std::string("H1"), H1FESpace_, true);

  gridfunctions_ = new hephaestus::GridFunctions;
  gridfunctions_->Register(std::string("source"), J_, true);
  gridfunctions_->Register(std::string("V"), V_, true);

  current_solver_options_ = new hephaestus::InputParameters;
  current_solver_options_->SetParam("Tolerance", float(1.0e-9));
  current_solver_options_->SetParam("MaxIter", (unsigned int)1000);
  current_solver_options_->SetParam("PrintLevel", 1);

  sps_params_ = new hephaestus::InputParameters;
  sps_params_->SetParam("SourceName", std::string("source"));
  sps_params_->SetParam("PotentialName", std::string("V"));
  sps_params_->SetParam("HCurlFESpaceName", std::string("HCurl"));
  sps_params_->SetParam("H1FESpaceName", std::string("H1"));
  sps_params_->SetParam("SolverOptions", *current_solver_options_);
  sps_params_->SetParam("ConductivityCoefName",
                        std::string("magnetic_permeability"));

  coefs_ = new hephaestus::Coefficients;
  coefs_->scalars.Register("magnetic_permeability", &coef1_, false);

  sps_ = new hephaestus::ScalarPotentialSource(*sps_params_);
  sps_->Init(*gridfunctions_, *fespaces_, *bc_maps_, *coefs_);

  mfem::ParLinearForm dummy(HCurlFESpace_);
  sps_->Apply(&dummy);

  // Normalise the current through the wedges and use them as a reference
  double flux = calcFlux(J_, ref_face_);
  *J_ /= abs(flux);
  if (V_)
    *V_ /= abs(flux);

  delete fespaces_;
  delete gridfunctions_;
  delete current_solver_options_;
  delete sps_params_;
  delete coefs_;
  delete sps_;
}

void OpenCoilSolver::setRefFace(const int face) { ref_face_ = face; }

} // namespace hephaestus