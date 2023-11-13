#include "closed_coil.hpp"

namespace hephaestus {

// Pushes an element into a vector if the vector does not yet contain that same
// element
template <typename T> void pushIfUnique(std::vector<T> &vec, const T el) {

  bool verify = true;

  for (auto e : vec) {
    if (e == el)
      verify = false;
  }

  if (verify == true)
    vec.push_back(el);
}

// Deletes and clears a vector of pointers
template <typename T> void deleteAndClear(std::vector<T *> v) {

  for (auto p : v)
    delete p;
  v.clear();
}

// Base class methods

ClosedCoilSolver::ClosedCoilSolver(const hephaestus::InputParameters &params,
                                   const mfem::Array<int> &coil_dom,
                                   const int electrode_face)
    : hcurl_fespace_name_(params.GetParam<std::string>("HCurlFESpaceName")),
      J_gf_name_(params.GetParam<std::string>("JGridFunctionName")),
      I_coef_name_(params.GetParam<std::string>("IFuncCoefName")),
      coil_domains_(coil_dom), mesh_parent_(nullptr), J_parent_(nullptr),
      HCurlFESpace_parent_(nullptr), m1_(nullptr) {

  elec_attrs_.first = electrode_face;
}

ClosedCoilSolver::~ClosedCoilSolver() {

  delete mesh_coil_;
  delete HCurlFESpace_coil_;
  delete H1FESpace_coil_;
  delete J_coil_;
  delete Jt_coil_;
}

void ClosedCoilSolver::Init(hephaestus::GridFunctions &gridfunctions,
                            const hephaestus::FESpaces &fespaces,
                            hephaestus::BCMap &bc_map,
                            hephaestus::Coefficients &coefficients) {

  // Retrieving the parent FE space and mesh
  HCurlFESpace_parent_ = fespaces.Get(hcurl_fespace_name_);
  if (HCurlFESpace_parent_ == nullptr) {
    const std::string error_message = hcurl_fespace_name_ +
                                      " not found in fespaces when "
                                      "creating ClosedCoilSolver\n";
    mfem::mfem_error(error_message.c_str());
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

  Itotal_ = coefficients.scalars.Get(I_coef_name_);
  if (Itotal_ == nullptr) {
    std::cout << I_coef_name_ + " not found in coefficients when "
                                "creating ClosedCoilSolver. "
                                "Assuming unit current. ";
    Itotal_ = new mfem::ConstantCoefficient(1.0);
  }

  mesh_parent_ = HCurlFESpace_parent_->GetParMesh();
  order_hcurl_ = HCurlFESpace_parent_->FEColl()->GetOrder();
  order_h1_ = order_hcurl_;

  makeWedge();
  prepareCoilSubmesh();
  solveTransition();
  solveCoil();
  buildM1();
  normaliseCurrent();
  restoreAttributes();
}

void ClosedCoilSolver::Apply(mfem::ParLinearForm *lf) {

  // The transformation and integration points themselves are not relevant, it's
  // just so we can call Eval
  mfem::ElementTransformation *Tr = mesh_parent_->GetElementTransformation(0);
  const mfem::IntegrationPoint &ip =
      mfem::IntRules.Get(J_coil_->ParFESpace()->GetFE(0)->GetGeomType(), 1)
          .IntPoint(0);

  double I = Itotal_->Eval(*Tr, ip);
  *J_coil_ *= I;
  mesh_coil_->Transfer(*J_coil_, *J_parent_);
  *J_coil_ /= I;

  m1_->Update();
  m1_->Assemble();
  m1_->AddMult(*J_parent_, *lf, 1.0);

}

void ClosedCoilSolver::SubtractSource(mfem::ParGridFunction *gf) {}

// ClosedCoilSolver main methods

void ClosedCoilSolver::makeWedge() {

  std::vector<int> bdr_els;

  // First we save the current domain attributes so they may be restored later
  for (int e = 0; e < mesh_parent_->GetNE(); ++e)
    old_dom_attrs.push_back(mesh_parent_->GetAttribute(e));

  new_domain_attr_ = mesh_parent_->attributes.Max() + 1;

  elec_attrs_.second = mesh_parent_->bdr_attributes.Max() + 1;

  // Now we need to find the electrode boundary
  for (int i = 0; i < mesh_parent_->GetNBE(); ++i) {
    if (mesh_parent_->GetBdrAttribute(i) == elec_attrs_.first) {
      bdr_els.push_back(i);
    }
  }

  Plane3D plane;

  if (bdr_els.size() > 0) {
    plane.make3DPlane(mesh_parent_, mesh_parent_->GetBdrFace(bdr_els[0]));
  }

  std::vector<int> elec_vtx;
  // Create a vector containing all of the vertices on the electrode
  for (auto b_fc : bdr_els) {

    mfem::Array<int> face_vtx;
    mesh_parent_->GetFaceVertices(mesh_parent_->GetBdrFace(b_fc), face_vtx);

    for (auto v : face_vtx)
      pushIfUnique(elec_vtx, v);
  }

  // Now we need to find all elements in the mesh that touch, on at least one
  // vertex, the electrode face if they do touch the vertex, are on one side of
  // the electrode, and belong to the coil domain, we add them to our wedge

  std::vector<int> wedge_els;

  for (int e = 0; e < mesh_parent_->GetNE(); ++e) {

    if (!isInDomain(e, coil_domains_, mesh_parent_) ||
        plane.side(elementCentre(e, mesh_parent_)) == 1)
      continue;

    mfem::Array<int> elem_vtx;
    mesh_parent_->GetElementVertices(e, elem_vtx);

    for (auto v1 : elem_vtx) {
      for (auto v2 : elec_vtx) {
        if (v1 == v2) {
          pushIfUnique(wedge_els, e);
        }
      }
    }
  }

  // Now we set the second electrode boundary attribute. Start with a list of
  // all the faces of the wedge elements and eliminate mesh and coil boundaries,
  // the first electrode, and faces between wedge elements

  std::vector<int> wedge_faces;
  mfem::Array<int> el_faces;
  mfem::Array<int> ori;

  for (auto e : wedge_els) {
    mesh_parent_->GetElementFaces(e, el_faces, ori);
    for (auto f : el_faces)
      pushIfUnique(wedge_faces, f);
  }

  for (auto wf : wedge_faces) {

    int e1, e2;
    mesh_parent_->GetFaceElements(wf, &e1, &e2);

    // If the face is a coil boundary
    if (!(isInDomain(e1, coil_domains_, mesh_parent_) &&
          isInDomain(e2, coil_domains_, mesh_parent_))) {
      continue;
    }

    // If the face is not true interior
    if (!(mesh_parent_->FaceIsInterior(wf) ||
          (mesh_parent_->GetFaceInformation(wf).tag ==
               mfem::Mesh::FaceInfoTag::SharedConforming ||
           mesh_parent_->GetFaceInformation(wf).tag ==
               mfem::Mesh::FaceInfoTag::SharedSlaveNonconforming))) {
      continue;
    }

    // If the face is shared between two elements internal to the wedge
    bool test1 = false;
    bool test2 = false;
    for (auto e : wedge_els) {
      if (e == e1)
        test1 = true;
      if (e == e2)
        test2 = true;
    }

    if (test1 && test2)
      continue;

    // If the face is part of the first electrode
    test1 = false;
    for (auto b_fc : bdr_els) {
      if (wf == mesh_parent_->GetBdrFace(b_fc)) {
        test1 = true;
        break;
      }
    }
    if (test1)
      continue;

    // At last, if the face is none of these things, it must be our second
    // electrode
    auto *new_elem = mesh_parent_->GetFace(wf)->Duplicate(mesh_parent_);
    new_elem->SetAttribute(elec_attrs_.second);
    mesh_parent_->AddBdrElement(new_elem);
  }

  // Only after this do we set the domain attributes
  for (auto e : wedge_els)
    mesh_parent_->SetAttribute(e, new_domain_attr_);

  transition_domain_.Append(new_domain_attr_);
  coil_domains_.Append(new_domain_attr_);

  mesh_parent_->FinalizeTopology();
  mesh_parent_->Finalize();
  mesh_parent_->SetAttributes();
}

void ClosedCoilSolver::prepareCoilSubmesh() {

  mesh_coil_ = new mfem::ParSubMesh(
      mfem::ParSubMesh::CreateFromDomain(*mesh_parent_, coil_domains_));

  inheritBdrAttributes(mesh_parent_, mesh_coil_);
  HCurlFESpace_coil_ = new mfem::ParFiniteElementSpace(
      mesh_coil_,
      new mfem::ND_FECollection(order_hcurl_, mesh_coil_->Dimension()));

  H1FESpace_coil_ = new mfem::ParFiniteElementSpace(
      mesh_coil_,
      new mfem::H1_FECollection(order_h1_, mesh_coil_->Dimension()));

  J_coil_ = new mfem::ParGridFunction(HCurlFESpace_coil_);
  Jt_coil_ = new mfem::ParGridFunction(HCurlFESpace_coil_);

  *J_coil_ = 0.0;
  *Jt_coil_ = 0.0;
}

void ClosedCoilSolver::solveTransition() {

  hephaestus::FESpaces fespaces;
  hephaestus::Coefficients coefs;
  hephaestus::BCMap bc_maps;

  hephaestus::GridFunctions gridfunctions;
  gridfunctions.Register("J_parent", J_parent_, false);

  hephaestus::InputParameters ocs_params;
  ocs_params.SetParam("SourceName", std::string("J_parent"));
  ocs_params.SetParam("IFuncCoefName", std::string("I"));
  ocs_params.SetParam("PotentialName", std::string("Phi"));

  hephaestus::OpenCoilSolver opencoil(ocs_params, transition_domain_,
                                      elec_attrs_);

  opencoil.Init(gridfunctions, fespaces, bc_maps, coefs);
  mfem::ParLinearForm dummy(HCurlFESpace_coil_);
  opencoil.Apply(&dummy);

  // The transition region result goes
  // Child -> Grandparent -> Parent
  // Ideally, it should go Child -> Parent
  // However, MFEM has issues creating transfer maps
  // between several generations
  mesh_coil_->Transfer(*J_parent_, *Jt_coil_);
}

void ClosedCoilSolver::solveCoil() {
  // (∇Va,∇ψ) = (Jt,∇ψ)
  // where Va is auxV_coil_, the auxiliary continuous "potential"
  // ψ are the H1 test functions
  // Jt is the vector function corresponding to the current in
  // the transition region
  // The boundary terms are zero because ∇Va and Jt are perpendicular
  // to the coil boundaries

  mfem::ParGridFunction auxV_coil(H1FESpace_coil_);
  auxV_coil = 0.0;

  mfem::ParBilinearForm a(H1FESpace_coil_);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  a.Assemble();

  mfem::VectorGridFunctionCoefficient JtCoef(Jt_coil_);
  mfem::ParLinearForm b(H1FESpace_coil_);
  b.AddDomainIntegrator(new mfem::DomainLFGradIntegrator(JtCoef));
  b.Assemble();

  mfem::HypreParMatrix A;
  mfem::Vector B, X;
  mfem::Array<int> boundary_dofs;
  a.FormLinearSystem(boundary_dofs, auxV_coil, b, A, X, B);

  mfem::HypreBoomerAMG amg(A);
  mfem::HyprePCG pcg(A);
  pcg.SetTol(1e-8);
  pcg.SetMaxIter(500);
  pcg.SetPrintLevel(1);
  pcg.SetPreconditioner(amg);
  pcg.Mult(B, X);
  a.RecoverFEMSolution(X, b, auxV_coil);

  // Now we form the final coil current
  mfem::ParDiscreteLinearOperator grad(H1FESpace_coil_, HCurlFESpace_coil_);
  grad.AddDomainInterpolator(new mfem::GradientInterpolator());
  grad.Assemble();

  grad.Mult(auxV_coil, *J_coil_);
  J_coil_->Add(-1.0, *Jt_coil_);
}

void ClosedCoilSolver::buildM1() {

  if (m1_ == nullptr)
    m1_ = new mfem::ParBilinearForm(HCurlFESpace_parent_);

  m1_->AddDomainIntegrator(
      new mfem::VectorFEMassIntegrator(new mfem::ConstantCoefficient(1.0)));
  m1_->Assemble();
}

void ClosedCoilSolver::normaliseCurrent() {

  double flux = calcFlux(J_coil_, elec_attrs_.first);
  *J_coil_ /= abs(flux);
}

void ClosedCoilSolver::restoreAttributes() {

  // Domain attributes
  for (int e = 0; e < mesh_parent_->GetNE(); ++e) {
    mesh_parent_->SetAttribute(e, old_dom_attrs[e]);
  }

  mesh_parent_->FinalizeTopology();
  mesh_parent_->Finalize();
  mesh_parent_->SetAttributes();
}

// Auxiliary methods

bool ClosedCoilSolver::isInDomain(const int el, const mfem::Array<int> &dom,
                                  const mfem::ParMesh *mesh) {

  // This is for ghost elements
  if (el < 0)
    return false;

  bool verify = false;

  for (auto sd : dom) {
    if (mesh->GetAttribute(el) == sd)
      verify = true;
  }

  return verify;
}

bool ClosedCoilSolver::isInDomain(const int el, const int &sd,
                                  const mfem::ParMesh *mesh) {

  // This is for ghost elements
  if (el < 0)
    return false;

  return mesh->GetAttribute(el) == sd;
}

mfem::Vector ClosedCoilSolver::elementCentre(int el, mfem::ParMesh *pm) {

  mfem::Array<int> elem_vtx;
  mfem::Vector com(3);
  com = 0.0;

  pm->GetElementVertices(el, elem_vtx);

  for (auto vtx : elem_vtx) {
    for (int j = 0; j < 3; ++j)
      com[j] += pm->GetVertex(vtx)[j] / (double)elem_vtx.Size();
  }

  return com;
}

// 3D Plane constructor and methods

Plane3D::Plane3D() : d(0) {

  u = new mfem::Vector(3);
  *u = 0.0;
}

Plane3D::~Plane3D() { delete u; }

void Plane3D::make3DPlane(const mfem::ParMesh *pm, const int face) {

  MFEM_ASSERT(pm->Dimension() == 3,
              "Plane3D only works in 3-dimensional meshes!");

  mfem::Array<int> face_vtx;
  std::vector<mfem::Vector> v;
  pm->GetFaceVertices(face, face_vtx);

  // First we get the coordinates of 3 vertices on the face
  for (auto vtx : face_vtx) {
    mfem::Vector vtx_coords(3);
    for (int j = 0; j < 3; ++j)
      vtx_coords[j] = pm->GetVertex(vtx)[j];
    v.push_back(vtx_coords);
  }

  // Now we find the unit vector normal to the face
  v[0] -= v[1];
  v[1] -= v[2];
  v[0].cross3D(v[1], *u);
  *u /= u->Norml2();

  // Finally, we find d:
  d = *u * v[2];
}

int Plane3D::side(const mfem::Vector v) {
  double val = *u * v - d;

  if (val > 0)
    return 1;
  else if (val < 0)
    return -1;
  else
    return 0;
}

}; // namespace hephaestus