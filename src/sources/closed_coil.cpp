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
                                   const int electrode_face, const int order)
    : hcurl_fespace_name_(params.GetParam<std::string>("HCurlFESpaceName")),
      J_gf_name_(params.GetParam<std::string>("JGridFunctionName")),
      I_coef_name_(params.GetParam<std::string>("IFuncCoefName")),
      coil_domains_(coil_dom), order_(order), coef1_(nullptr), coef0_(nullptr),
      mesh_parent_(nullptr), J_parent_(nullptr), HCurlFESpace_parent_(nullptr) {

  elec_attrs_.first = electrode_face;
  coef1_ = new mfem::ConstantCoefficient(1.0);
  coef0_ = new mfem::ConstantCoefficient(0.0);
}

ClosedCoilSolver::~ClosedCoilSolver() {

  delete coef1_;
  delete coef0_;

  // Deal with destructor!!
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
                                      "creating ClosedCoilSolver\n";
    mfem::mfem_error(error_message.c_str());
  }

  Itotal_ = coefficients.scalars.Get(I_coef_name_);
  if (Itotal_ == nullptr) {
    std::cout << I_coef_name_ + " not found in coefficients when "
                                "creating ClosedCoilSolver. "
                                "Assuming unit current. ";
    Itotal_ = new mfem::ConstantCoefficient(1.0);
  }

  mesh_parent_ = HCurlFESpace_parent_->GetParMesh();

  makeWedge();
  prepareCoilSubmesh();
  solveTransition();
  restoreAttributes();
}

void ClosedCoilSolver::Apply(mfem::ParLinearForm *lf) {

  mesh_coil_->Transfer(*J_coil_, *J_parent_);
  lf->Add(1.0, *J_parent_);
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

void ClosedCoilSolver::cleanDivergence(hephaestus::GridFunctions *gridfunctions,
                                       std::string J_name, std::string V_name,
                                       hephaestus::BCMap *bc_map) {

  hephaestus::InputParameters pars;
  hephaestus::FESpaces fes;

  pars.SetParam("VectorGridFunctionName", J_name);
  pars.SetParam("ScalarGridFunctionName", V_name);
  hephaestus::HelmholtzProjector projector(pars);
  projector.Project(*gridfunctions, fes, *bc_map);
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

void ClosedCoilSolver::solveTransition() {

  ocs_params_ = new hephaestus::InputParameters;
  bc_maps_ = new hephaestus::BCMap;
  coefs_ = new hephaestus::Coefficients;
  fespaces_ = new hephaestus::FESpaces;
  gridfunctions_ = new hephaestus::GridFunctions;
  gridfunctions_->Register("J_coil", J_coil_, false);

  ocs_params_->SetParam("SourceName", std::string("J_coil"));
  ocs_params_->SetParam("IFuncCoefName", std::string("I"));
  ocs_params_->SetParam("PotentialName", std::string("Phi"));

  opencoil_ = new hephaestus::OpenCoilSolver(*ocs_params_, transition_domain_,
                                             elec_attrs_, order_);
  opencoil_->Init(*gridfunctions_, *fespaces_, *bc_maps_, *coefs_);
  mfem::ParLinearForm dummy;
  opencoil_->Apply(&dummy);
}

void ClosedCoilSolver::prepareCoilSubmesh() {

  // Extracting submesh
  mesh_coil_ = new mfem::ParSubMesh(
      mfem::ParSubMesh::CreateFromDomain(*mesh_parent_, coil_domains_));

  inheritBdrAttributes(mesh_parent_, mesh_coil_);
  
  // FES and GridFunction
  HCurlFESpace_coil_ = new mfem::ParFiniteElementSpace(
      mesh_coil_, new mfem::ND_FECollection(order_, mesh_coil_->Dimension()));
  J_coil_ = new mfem::ParGridFunction(HCurlFESpace_coil_);
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

void ClosedCoilSolver::subdomainToArray(
    const std::vector<hephaestus::Subdomain> &sd, mfem::Array<int> &arr) {

  arr.DeleteAll();
  for (auto s : sd)
    arr.Append(s.id);
}

void ClosedCoilSolver::subdomainToArray(const hephaestus::Subdomain &sd,
                                        mfem::Array<int> &arr) {

  arr.DeleteAll();
  arr.Append(sd.id);
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