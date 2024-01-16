#include "closed_coil.hpp"

namespace hephaestus
{

// Pushes an element into a vector if the vector does not yet contain that same
// element
template <typename T>
void
pushIfUnique(std::vector<T> & vec, const T el)
{

  bool verify = true;

  for (auto e : vec)
  {
    if (e == el)
      verify = false;
  }

  if (verify == true)
    vec.push_back(el);
}

// Deletes and clears a vector of pointers
template <typename T>
void
deleteAndClear(std::vector<T *> v)
{

  for (auto p : v)
    delete p;
  v.clear();
}

// Base class methods

ClosedCoilSolver::ClosedCoilSolver(const hephaestus::InputParameters & params,
                                   const mfem::Array<int> & coil_dom,
                                   const int electrode_face)
  : hcurl_fespace_name_(params.GetParam<std::string>("HCurlFESpaceName")),
    h1_fespace_name_(params.GetParam<std::string>("H1FESpaceName")),
    J_gf_name_(params.GetParam<std::string>("JGridFunctionName")),
    I_coef_name_(params.GetParam<std::string>("IFuncCoefName")),
    J_transfer_(params.GetOptionalParam<bool>("JTransfer", false)),
    coil_domains_(coil_dom)
{

  hephaestus::InputParameters default_pars;
  default_pars.SetParam("Tolerance", float(1e-18));
  default_pars.SetParam("AbsTolerance", float(1e-18));
  default_pars.SetParam("MaxIter", (unsigned int)1000);
  default_pars.SetParam("PrintLevel", 1);

  solver_options_ =
      params.GetOptionalParam<hephaestus::InputParameters>("SolverOptions", default_pars);

  elec_attrs_.first = electrode_face;
}

ClosedCoilSolver::~ClosedCoilSolver()
{
  if (owns_Itotal_)
    delete Itotal_;
}

void
ClosedCoilSolver::Init(hephaestus::GridFunctions & gridfunctions,
                       const hephaestus::FESpaces & fespaces,
                       hephaestus::BCMap & bc_map,
                       hephaestus::Coefficients & coefficients)
{

  // Retrieving the parent FE space and mesh

  HCurlFESpace_parent_ = fespaces.Get(hcurl_fespace_name_);
  if (HCurlFESpace_parent_ == nullptr)
  {
    const std::string error_message = hcurl_fespace_name_ + " not found in fespaces when "
                                                            "creating ClosedCoilSolver\n";
    mfem::mfem_error(error_message.c_str());
  }

  mesh_parent_ = HCurlFESpace_parent_->GetParMesh();
  order_hcurl_ = HCurlFESpace_parent_->FEColl()->GetOrder();
  order_h1_ = order_hcurl_;

  // Optional FE Spaces and parameters
  H1FESpace_parent_ = fespaces.Get(h1_fespace_name_);

  if (H1FESpace_parent_ == nullptr)
  {
    std::cout << h1_fespace_name_ + " not found in fespaces when "
                                    "creating ClosedCoilSolver. Creating from mesh.\n";

    // Need to free this memory after use. FEC not freed by
    // ParFiniteElementSpace destructor!
    H1FESpace_parent_fec_ =
        std::make_unique<mfem::H1_FECollection>(order_h1_, mesh_parent_->Dimension());

    H1FESpace_parent_ = new mfem::ParFiniteElementSpace(mesh_parent_, H1FESpace_parent_fec_.get());
  }

  J_parent_ = gridfunctions.Get(J_gf_name_);
  if (J_parent_ == nullptr)
  {
    std::cout << J_gf_name_ + " not found in gridfunctions when "
                              "creating OpenCoilSolver. Creating new GridFunction.\n";

    J_parent_ = new mfem::ParGridFunction(HCurlFESpace_parent_);
  }
  else if (J_parent_->ParFESpace()->FEColl()->GetContType() !=
           mfem::FiniteElementCollection::TANGENTIAL)
  {
    mfem::mfem_error("J GridFunction must be of HCurl type.");
  }

  Itotal_ = coefficients.scalars.Get(I_coef_name_);
  if (Itotal_ == nullptr)
  {
    std::cout << I_coef_name_ + " not found in coefficients when "
                                "creating ClosedCoilSolver. "
                                "Assuming unit current. ";

    Itotal_ = new mfem::ConstantCoefficient(1.0);
    owns_Itotal_ = true; // Responsible for deleting this.
  }

  if (final_lf_ == nullptr)
  {
    final_lf_ = std::make_unique<mfem::ParLinearForm>(HCurlFESpace_parent_);
    *final_lf_ = 0.0;
  }

  makeWedge();
  prepareCoilSubmesh();
  solveTransition();
  solveCoil();
  restoreAttributes();

  // Cleanup memory.
  if (!fespaces.Has(h1_fespace_name_))
  {
    delete H1FESpace_parent_;
    H1FESpace_parent_fec_.reset();
  }

  if (!gridfunctions.Has(J_gf_name_))
    delete J_parent_;
}

void
ClosedCoilSolver::Apply(mfem::ParLinearForm * lf)
{

  // The transformation and integration points themselves are not relevant, it's
  // just so we can call Eval
  mfem::ElementTransformation * Tr = mesh_parent_->GetElementTransformation(0);
  const mfem::IntegrationPoint & ip =
      mfem::IntRules.Get(HCurlFESpace_parent_->GetFE(0)->GetGeomType(), 1).IntPoint(0);

  double I = Itotal_->Eval(*Tr, ip);
  lf->Add(I, *final_lf_);

  if (J_transfer_)
  {
    *J_parent_ = 0.0;
    J_parent_->Add(I, *Jt_parent_);
  }
}

void
ClosedCoilSolver::SubtractSource(mfem::ParGridFunction * gf)
{
}

// ClosedCoilSolver main methods

void
ClosedCoilSolver::makeWedge()
{

  std::vector<int> bdr_els;

  // First we save the current domain attributes so they may be restored later
  for (int e = 0; e < mesh_parent_->GetNE(); ++e)
    old_dom_attrs.push_back(mesh_parent_->GetAttribute(e));

  new_domain_attr_ = mesh_parent_->attributes.Max() + 1;

  elec_attrs_.second = mesh_parent_->bdr_attributes.Max() + 1;

  // Now we need to find the electrode boundary
  for (int i = 0; i < mesh_parent_->GetNBE(); ++i)
  {
    if (mesh_parent_->GetBdrAttribute(i) == elec_attrs_.first)
    {
      bdr_els.push_back(i);
    }
  }

  Plane3D plane;

  if (bdr_els.size() > 0)
  {
    plane.make3DPlane(mesh_parent_, mesh_parent_->GetBdrElementFaceIndex(bdr_els[0]));
  }

  std::vector<int> elec_vtx;
  // Create a vector containing all of the vertices on the electrode
  for (auto b_fc : bdr_els)
  {

    mfem::Array<int> face_vtx;
    mesh_parent_->GetFaceVertices(mesh_parent_->GetBdrElementFaceIndex(b_fc), face_vtx);

    for (auto v : face_vtx)
      pushIfUnique(elec_vtx, v);
  }

  // Now we need to find all elements in the mesh that touch, on at least one
  // vertex, the electrode face if they do touch the vertex, are on one side of
  // the electrode, and belong to the coil domain, we add them to our wedge

  std::vector<int> wedge_els;

  for (int e = 0; e < mesh_parent_->GetNE(); ++e)
  {

    if (!isInDomain(e, coil_domains_, mesh_parent_) ||
        plane.side(elementCentre(e, mesh_parent_)) == 1)
      continue;

    mfem::Array<int> elem_vtx;
    mesh_parent_->GetElementVertices(e, elem_vtx);

    for (auto v1 : elem_vtx)
    {
      for (auto v2 : elec_vtx)
      {
        if (v1 == v2)
        {
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

  for (auto e : wedge_els)
  {
    mesh_parent_->GetElementFaces(e, el_faces, ori);
    for (auto f : el_faces)
      pushIfUnique(wedge_faces, f);
  }

  for (auto wf : wedge_faces)
  {

    int e1, e2;
    mesh_parent_->GetFaceElements(wf, &e1, &e2);

    // If the face is a coil boundary
    if (!(isInDomain(e1, coil_domains_, mesh_parent_) &&
          isInDomain(e2, coil_domains_, mesh_parent_)))
    {
      continue;
    }

    // If the face is not true interior
    if (!(mesh_parent_->FaceIsInterior(wf) ||
          (mesh_parent_->GetFaceInformation(wf).tag == mfem::Mesh::FaceInfoTag::SharedConforming ||
           mesh_parent_->GetFaceInformation(wf).tag ==
               mfem::Mesh::FaceInfoTag::SharedSlaveNonconforming)))
    {
      continue;
    }

    // If the face is shared between two elements internal to the wedge
    bool test1 = false;
    bool test2 = false;
    for (auto e : wedge_els)
    {
      if (e == e1)
        test1 = true;
      if (e == e2)
        test2 = true;
    }

    if (test1 && test2)
      continue;

    // If the face is part of the first electrode
    test1 = false;
    for (auto b_fc : bdr_els)
    {
      if (wf == mesh_parent_->GetBdrElementFaceIndex(b_fc))
      {
        test1 = true;
        break;
      }
    }
    if (test1)
      continue;

    // At last, if the face is none of these things, it must be our second
    // electrode
    auto * new_elem = mesh_parent_->GetFace(wf)->Duplicate(mesh_parent_);
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

void
ClosedCoilSolver::prepareCoilSubmesh()
{

  mesh_coil_ = std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromDomain(*mesh_parent_, coil_domains_));

  H1FESpace_coil_fec_ = std::make_unique<mfem::H1_FECollection>(order_h1_, mesh_coil_->Dimension());

  H1FESpace_coil_ =
      std::make_unique<mfem::ParFiniteElementSpace>(mesh_coil_.get(), H1FESpace_coil_fec_.get());

  Jaux_coil_fec_ = std::make_unique<mfem::ND_FECollection>(order_hcurl_, mesh_coil_->Dimension());

  Jaux_coil_ = std::make_unique<mfem::ParGridFunction>(
      new mfem::ParFiniteElementSpace(mesh_coil_.get(), Jaux_coil_fec_.get()));
  *Jaux_coil_ = 0.0;

  V_coil_ = std::make_unique<mfem::ParGridFunction>(H1FESpace_coil_.get());
  *V_coil_ = 0.0;

  mesh_t_ = std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromDomain(*mesh_parent_, transition_domain_));
}

void
ClosedCoilSolver::solveTransition()
{

  mfem::ParGridFunction V_parent(H1FESpace_parent_);
  V_parent = 0.0;

  hephaestus::FESpaces fespaces;
  hephaestus::Coefficients coefs;
  hephaestus::BCMap bc_maps;

  hephaestus::GridFunctions gridfunctions;
  gridfunctions.Register("J_parent", J_parent_, false);
  gridfunctions.Register("V_parent", &V_parent, false);

  hephaestus::InputParameters ocs_params;
  ocs_params.SetParam("SourceName", std::string("J_parent"));
  ocs_params.SetParam("IFuncCoefName", std::string("I"));
  ocs_params.SetParam("PotentialName", std::string("V_parent"));
  ocs_params.SetParam("SolverOptions", solver_options_);

  hephaestus::OpenCoilSolver opencoil(ocs_params, transition_domain_, elec_attrs_);

  opencoil.Init(gridfunctions, fespaces, bc_maps, coefs);
  opencoil.Apply(final_lf_.get());

  mesh_coil_->Transfer(V_parent, *V_coil_);
}

void
ClosedCoilSolver::solveCoil()
{
  // (∇Va,∇ψ) = (∇Vt,∇ψ)
  // where Va is Vaux_coil_, the auxiliary continuous "potential"
  // ψ are the H1 test functions
  // Vt is the transition potential
  // The boundary terms are zero because ∇Va and ∇Vt are perpendicular
  // to the coil boundaries

  mfem::ParGridFunction Vaux_coil(H1FESpace_coil_.get());
  Vaux_coil = 0.0;

  mfem::ParBilinearForm a_t(H1FESpace_coil_.get());
  mfem::ParLinearForm b_coil(H1FESpace_coil_.get());
  b_coil = 0.0;

  attrToMarker(transition_domain_, transition_markers_, mesh_coil_->attributes.Max());
  a_t.AddDomainIntegrator(new mfem::DiffusionIntegrator, transition_markers_);
  a_t.Assemble();
  a_t.Finalize();
  a_t.AddMult(*V_coil_, b_coil, 1.0);

  mfem::ParBilinearForm a_coil(H1FESpace_coil_.get());
  a_coil.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  a_coil.Assemble();

  mfem::Array<int> ess_bdr_tdofs_coil;

  // This creates a binary representation of which MPI ranks contain at
  // least one element
  int ref_rank = 0;
  int has_els = (bool)mesh_coil_->GetNE() ? 1 << mfem::Mpi::WorldRank() : 0;
  int has_els_sum;

  MPI_Allreduce(&has_els, &has_els_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MFEM_ASSERT(has_els_sum != 0, "Empty coil submesh!");

  for (int i = 0; i < mfem::Mpi::WorldSize(); ++i)
  {
    if ((1 << i & has_els_sum) != 0)
    {
      ref_rank = i;
      break;
    }
  }

  if (mfem::Mpi::WorldRank() == ref_rank)
  {
    ess_bdr_tdofs_coil.SetSize(1);
    ess_bdr_tdofs_coil[0] = 0;
  }

  mfem::HypreParMatrix A0_coil;
  mfem::Vector X0_coil;
  mfem::Vector B0_coil;
  a_coil.FormLinearSystem(ess_bdr_tdofs_coil, Vaux_coil, b_coil, A0_coil, X0_coil, B0_coil);
  hephaestus::DefaultH1PCGSolver a_coil_solver(solver_options_, A0_coil);
  a_coil_solver.Mult(B0_coil, X0_coil);
  a_coil.RecoverFEMSolution(X0_coil, b_coil, Vaux_coil);

  // Now we form the final coil current
  mfem::ParDiscreteLinearOperator grad(H1FESpace_coil_.get(), Jaux_coil_->ParFESpace());
  grad.AddDomainInterpolator(new mfem::GradientInterpolator());
  grad.Assemble();
  grad.Mult(Vaux_coil, *Jaux_coil_);

  if (J_transfer_)
    Jt_parent_ = std::make_unique<mfem::ParGridFunction>(*J_parent_);

  *J_parent_ = 0.0;
  mesh_coil_->Transfer(*Jaux_coil_, *J_parent_);

  mfem::ParBilinearForm m1(HCurlFESpace_parent_);
  hephaestus::attrToMarker(coil_domains_, coil_markers_, mesh_parent_->attributes.Max());

  auto constant_coeff = std::make_unique<mfem::ConstantCoefficient>(1.0);
  m1.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(constant_coeff.get()), coil_markers_);
  m1.Assemble();
  m1.AddMult(*J_parent_, *final_lf_, -1.0);

  // We can't properly calculate the flux of Jaux on the parent mesh, so we
  // transfer it first to the transition mesh. This will be used in the
  // normalisation step
  Jaux_t_fec_ = std::make_unique<mfem::ND_FECollection>(
      mfem::ND_FECollection(order_hcurl_, mesh_t_->Dimension()));

  mfem::ParGridFunction Jaux_t_(new mfem::ParFiniteElementSpace(mesh_t_.get(), Jaux_t_fec_.get()));
  Jaux_t_ = 0.0;

  mesh_t_->Transfer(*J_parent_, Jaux_t_);

  // The total flux across the electrode face is Φ_t-Φ_aux
  // where Φ_t is the transition flux, already normalised to be 1
  double flux = 1.0 - calcFlux(&Jaux_t_, elec_attrs_.first);

  if (J_transfer_)
  {
    *J_parent_ -= *Jt_parent_;
    *J_parent_ /= -flux;
    *Jt_parent_ = *J_parent_;
  }

  *final_lf_ /= flux;
}

void
ClosedCoilSolver::restoreAttributes()
{

  // Domain attributes
  for (int e = 0; e < mesh_parent_->GetNE(); ++e)
  {
    mesh_parent_->SetAttribute(e, old_dom_attrs[e]);
  }

  mesh_parent_->FinalizeTopology();
  mesh_parent_->Finalize();
  mesh_parent_->SetAttributes();
}

// Auxiliary methods

bool
ClosedCoilSolver::isInDomain(const int el, const mfem::Array<int> & dom, const mfem::ParMesh * mesh)
{

  // This is for ghost elements
  if (el < 0)
    return false;

  bool verify = false;

  for (auto sd : dom)
  {
    if (mesh->GetAttribute(el) == sd)
      verify = true;
  }

  return verify;
}

bool
ClosedCoilSolver::isInDomain(const int el, const int & sd, const mfem::ParMesh * mesh)
{

  // This is for ghost elements
  if (el < 0)
    return false;

  return mesh->GetAttribute(el) == sd;
}

mfem::Vector
ClosedCoilSolver::elementCentre(int el, mfem::ParMesh * pm)
{

  mfem::Array<int> elem_vtx;
  mfem::Vector com(3);
  com = 0.0;

  pm->GetElementVertices(el, elem_vtx);

  for (auto vtx : elem_vtx)
  {
    for (int j = 0; j < 3; ++j)
      com[j] += pm->GetVertex(vtx)[j] / (double)elem_vtx.Size();
  }

  return com;
}

// 3D Plane constructor and methods

Plane3D::Plane3D() : d(0)
{

  u = std::make_unique<mfem::Vector>(3);
  *u = 0.0;
}

void
Plane3D::make3DPlane(const mfem::ParMesh * pm, const int face)
{

  MFEM_ASSERT(pm->Dimension() == 3, "Plane3D only works in 3-dimensional meshes!");

  mfem::Array<int> face_vtx;
  std::vector<mfem::Vector> v;
  pm->GetFaceVertices(face, face_vtx);

  // First we get the coordinates of 3 vertices on the face
  for (auto vtx : face_vtx)
  {
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

int
Plane3D::side(const mfem::Vector v)
{
  double val = *u * v - d;

  if (val > 0)
    return 1;
  else if (val < 0)
    return -1;
  else
    return 0;
}

}; // namespace hephaestus