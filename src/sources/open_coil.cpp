#include "open_coil.hpp"

namespace hephaestus {

///// THESE FUNCTIONS WILL EVENTUALLY GO INTO A UTILS FILE ///////////

double highV(const mfem::Vector &x, double t) { return 0.5; }
double lowV(const mfem::Vector &x, double t) { return -0.5; }

double calcFlux(mfem::GridFunction *v_field, int face_attr,
                mfem::Coefficient &q) {

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
        mesh->GetFaceElementTransformations(mesh->GetBdrElementFaceIndex(i));
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
      flux += q.Eval(*FTr, ip) * val * ip.weight * face_weight;
    }
  }

  double total_flux;
  MPI_Allreduce(&flux, &total_flux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return total_flux;
}

double calcFlux(mfem::GridFunction *v_field, int face_attr) {
  mfem::ConstantCoefficient one_coef(1.0);
  return calcFlux(v_field, face_attr, one_coef);
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

void attrToMarker(const mfem::Array<int> attr_list,
                  mfem::Array<int> &marker_list, int max_attr) {

  marker_list.SetSize(max_attr);
  marker_list = 0;

  for (auto a : attr_list)
    marker_list[a - 1] = 1;
}

void cleanDivergence(mfem::ParGridFunction &Vec_GF,
                     hephaestus::InputParameters solve_pars) {

  hephaestus::InputParameters pars;
  hephaestus::GridFunctions gfs;
  hephaestus::FESpaces fes;
  hephaestus::BCMap bcs;

  gfs.Register("Vector_GF", &Vec_GF, false);
  pars.SetParam("VectorGridFunctionName", std::string("Vector_GF"));
  pars.SetParam("SolverOptions", solve_pars);
  hephaestus::HelmholtzProjector projector(pars);
  projector.Project(gfs, fes, bcs);
}

void cleanDivergence(hephaestus::GridFunctions &gfs, hephaestus::BCMap &bcs,
                     const std::string vec_gf_name,
                     const std::string scalar_gf_name,
                     hephaestus::InputParameters solve_pars) {

  hephaestus::InputParameters pars;
  hephaestus::FESpaces fes;

  pars.SetParam("VectorGridFunctionName", vec_gf_name);
  pars.SetParam("ScalarGridFunctionName", scalar_gf_name);
  pars.SetParam("SolverOptions", solve_pars);
  hephaestus::HelmholtzProjector projector(pars);
  projector.Project(gfs, fes, bcs);
}

/////////////////////////////////////////////////////////////////////

OpenCoilSolver::OpenCoilSolver(const hephaestus::InputParameters &params,
                               const mfem::Array<int> &coil_dom,
                               const std::pair<int, int> electrodes)
    : grad_phi_name_(params.GetParam<std::string>("GradPotentialName")),
      V_gf_name_(params.GetParam<std::string>("PotentialName")),
      I_coef_name_(params.GetParam<std::string>("IFuncCoefName")),
      cond_coef_name_(params.GetParam<std::string>("ConductivityCoefName")),
      coil_domains_(coil_dom), elec_attrs_(electrodes), sigma_(nullptr),
      mesh_parent_(nullptr), mesh_(nullptr), H1FESpace_(nullptr),
      HCurlFESpace_(nullptr), grad_p_parent_(nullptr), grad_p_t_parent_(nullptr),
      V_parent_(nullptr), Vt_parent_(nullptr), grad_p_(nullptr), V_(nullptr),
      m1_(nullptr), final_lf_(nullptr), high_src_(highV), low_src_(lowV),
      high_terminal_(1), low_terminal_(1), grad_phi_transfer_(true) {

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
  ifDelete(grad_p_);
  ifDelete(V_);
  ifDelete(grad_p_t_parent_);
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

  sigma_ = coefficients.scalars.Get(cond_coef_name_);
  if (sigma_ == nullptr) {
    std::cout << cond_coef_name_ + " not found in coefficients when "
                                   "creating OpenCoilSolver. "
                                   "Assuming unit conductivity.\n";
    std::cout << "Warning: GradPhi field undefined. The GridFunction "
                 "associated with it will be set to zero.\n";
    sigma_ = new mfem::ConstantCoefficient(1.0);
    grad_phi_transfer_ = false;
  }

  grad_p_parent_ = gridfunctions.Get(grad_phi_name_);
  if (grad_p_parent_ == nullptr) {
    const std::string error_message = grad_phi_name_ +
                                      " not found in gridfunctions when "
                                      "creating OpenCoilSolver\n";
    mfem::mfem_error(error_message.c_str());
  } else if (grad_p_parent_->ParFESpace()->FEColl()->GetContType() !=
             mfem::FiniteElementCollection::TANGENTIAL) {
    mfem::mfem_error("J GridFunction must be of HCurl type.");
  }
  order_hcurl_ = grad_p_parent_->ParFESpace()->FEColl()->GetOrder();

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

  mesh_parent_ = grad_p_parent_->ParFESpace()->GetParMesh();

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
      mfem::IntRules.Get(grad_p_parent_->ParFESpace()->GetFE(0)->GetGeomType(), 1)
          .IntPoint(0);

  double I = Itotal_->Eval(*Tr, ip);

  *grad_p_parent_ = 0.0;
  if (grad_phi_transfer_)
    grad_p_parent_->Add(I, *grad_p_t_parent_);
  
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

  if (grad_p_ == nullptr)
    grad_p_ = new mfem::ParGridFunction(HCurlFESpace_);

  if (grad_p_t_parent_ == nullptr)
    grad_p_t_parent_ = new mfem::ParGridFunction(*grad_p_parent_);

  *V_ = 0.0;
  *grad_p_ = 0.0;
  *grad_p_t_parent_ = 0.0;
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
  fespaces.Register(std::string("HCurl"), HCurlFESpace_, false);
  fespaces.Register(std::string("H1"), H1FESpace_, false);

  hephaestus::GridFunctions gridfunctions;
  gridfunctions.Register(std::string("GradPhi"), grad_p_, false);
  gridfunctions.Register(std::string("V"), V_, false);

  hephaestus::InputParameters sps_params;
  sps_params.SetParam("GradPotentialName", std::string("GradPhi"));
  sps_params.SetParam("PotentialName", std::string("V"));
  sps_params.SetParam("HCurlFESpaceName", std::string("HCurl"));
  sps_params.SetParam("H1FESpaceName", std::string("H1"));
  sps_params.SetParam("SolverOptions", solver_options_);
  sps_params.SetParam("ConductivityCoefName",
                      std::string("electric_conductivity"));

  hephaestus::Coefficients coefs;
  coefs.scalars.Register("electric_conductivity", sigma_, false);

  hephaestus::ScalarPotentialSource sps(sps_params);
  sps.Init(gridfunctions, fespaces, bc_maps, coefs);

  mfem::ParLinearForm dummy(HCurlFESpace_);
  sps.Apply(&dummy);

  // Normalise the current through the wedges and use them as a reference
  double flux = calcFlux(grad_p_, ref_face_, *sigma_);
  *grad_p_ /= abs(flux);
  if (V_)
    *V_ /= abs(flux);

  mesh_->Transfer(*grad_p_, *grad_p_t_parent_);
  if (V_parent_)
    mesh_->Transfer(*V_, *Vt_parent_);

  buildM1();

  final_lf_ = new mfem::ParLinearForm(grad_p_t_parent_->ParFESpace());
  *final_lf_ = 0.0;
  m1_->AddMult(*grad_p_t_parent_, *final_lf_, 1.0);
}

void OpenCoilSolver::buildM1() {

  if (m1_ == nullptr) {

    m1_ = new mfem::ParBilinearForm(grad_p_parent_->ParFESpace());
    hephaestus::attrToMarker(coil_domains_, coil_markers_,
                             mesh_parent_->attributes.Max());
    m1_->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(sigma_),
                             coil_markers_);
    m1_->Assemble();
    m1_->Finalize();
  }
}

void OpenCoilSolver::setRefFace(const int face) { ref_face_ = face; }

} // namespace hephaestus