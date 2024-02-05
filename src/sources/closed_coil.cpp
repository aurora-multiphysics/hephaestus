#include "closed_coil.hpp"

#include <utility>

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

// Base class methods

ClosedCoilSolver::ClosedCoilSolver(const hephaestus::InputParameters & params,
                                   mfem::Array<int> coil_dom,
                                   const int electrode_face)
  : _hcurl_fespace_name(params.GetParam<std::string>("HCurlFESpaceName")),
    _h1_fespace_name(params.GetParam<std::string>("H1FESpaceName")),
    _grad_phi_name(params.GetParam<std::string>("GradPotentialName")),
    _i_coef_name(params.GetParam<std::string>("IFuncCoefName")),
    _cond_coef_name(params.GetParam<std::string>("ConductivityCoefName")),
    _grad_phi_transfer(params.GetOptionalParam<bool>("GradPhiTransfer", false)),
    _coil_domains(std::move(coil_dom))
{
  hephaestus::InputParameters default_pars;
  default_pars.SetParam("Tolerance", float(1e-18));
  default_pars.SetParam("AbsTolerance", float(1e-18));
  default_pars.SetParam("MaxIter", (unsigned int)1000);
  default_pars.SetParam("PrintLevel", 1);

  _solver_options =
      params.GetOptionalParam<hephaestus::InputParameters>("SolverOptions", default_pars);

  _elec_attrs.first = electrode_face;
}

void
ClosedCoilSolver::Init(hephaestus::GridFunctions & gridfunctions,
                       const hephaestus::FESpaces & fespaces,
                       hephaestus::BCMap & bc_map,
                       hephaestus::Coefficients & coefficients)
{

  // Retrieving the parent FE space and mesh

  _h_curl_fe_space_parent = fespaces.Get(_hcurl_fespace_name);

  _mesh_parent = _h_curl_fe_space_parent->GetParMesh();
  _order_hcurl = _h_curl_fe_space_parent->FEColl()->GetOrder();
  _order_h1 = _order_hcurl;

  // Optional FE Spaces and parameters
  _h1_fe_space_parent = fespaces.Get(_h1_fespace_name);
  if (_h1_fe_space_parent == nullptr)
  {
    std::cout << _h1_fespace_name + " not found in fespaces when "
                                    "creating ClosedCoilSolver. Creating from mesh.\n";

    // Need to free this memory after use. FEC not freed by ParFiniteElementSpace destructor!
    _h1_fe_space_parent_fec =
        std::make_unique<mfem::H1_FECollection>(_order_h1, _mesh_parent->Dimension());

    _h1_fe_space_parent =
        std::make_shared<mfem::ParFiniteElementSpace>(_mesh_parent, _h1_fe_space_parent_fec.get());
  }

  _grad_phi_parent = gridfunctions.Get(_grad_phi_name);
  if (_grad_phi_parent == nullptr)
  {
    std::cout << _grad_phi_name + " not found in gridfunctions when "
                                  "creating OpenCoilSolver. Creating new GridFunction.\n";
    _grad_phi_parent = std::make_shared<mfem::ParGridFunction>(_h_curl_fe_space_parent.get());
  }
  else if (_grad_phi_parent->ParFESpace()->FEColl()->GetContType() !=
           mfem::FiniteElementCollection::TANGENTIAL)
  {
    mfem::mfem_error("GradPhi GridFunction must be of HCurl type.");
  }

  _itotal = coefficients._scalars.Get(_i_coef_name);
  if (_itotal == nullptr)
  {
    std::cout << _i_coef_name + " not found in coefficients when "
                                "creating ClosedCoilSolver. "
                                "Assuming unit current. ";

    _itotal = std::make_shared<mfem::ConstantCoefficient>(1.0);
  }

  _sigma = coefficients._scalars.Get(_cond_coef_name);
  if (_sigma == nullptr)
  {
    std::cout << _cond_coef_name + " not found in coefficients when "
                                   "creating ClosedCoilSolver. "
                                   "Assuming unit conductivity.\n";
    std::cout << "Warning: GradPhi field undefined. The GridFunction "
                 "associated with it will be set to zero.\n";

    _sigma = std::make_shared<mfem::ConstantCoefficient>(1.0);

    _grad_phi_transfer = false;
  }

  if (_final_lf == nullptr)
  {
    _final_lf = std::make_unique<mfem::ParLinearForm>(_h_curl_fe_space_parent.get());
    *_final_lf = 0.0;
  }

  MakeWedge();
  PrepareCoilSubmesh();
  SolveTransition();
  SolveCoil();
  RestoreAttributes();
}

void
ClosedCoilSolver::Apply(mfem::ParLinearForm * lf)
{

  // The transformation and integration points themselves are not relevant, it's
  // just so we can call Eval
  mfem::ElementTransformation * tr = _mesh_parent->GetElementTransformation(0);
  const mfem::IntegrationPoint & ip =
      mfem::IntRules.Get(_h_curl_fe_space_parent->GetFE(0)->GetGeomType(), 1).IntPoint(0);

  double i = _itotal->Eval(*tr, ip);
  lf->Add(i, *_final_lf);

  *_grad_phi_parent = 0.0;
  if (_grad_phi_transfer)
  {
    _grad_phi_parent->Add(i, *_grad_phi_t_parent);
  }
}

void
ClosedCoilSolver::SubtractSource(mfem::ParGridFunction * gf)
{
}

// ClosedCoilSolver main methods

void
ClosedCoilSolver::MakeWedge()
{
  std::vector<int> bdr_els;

  // First we save the current domain attributes so they may be restored later
  for (int e = 0; e < _mesh_parent->GetNE(); ++e)
    _old_dom_attrs.push_back(_mesh_parent->GetAttribute(e));

  _new_domain_attr = _mesh_parent->attributes.Max() + 1;

  _elec_attrs.second = _mesh_parent->bdr_attributes.Max() + 1;

  // Now we need to find the electrode boundary
  for (int i = 0; i < _mesh_parent->GetNBE(); ++i)
  {
    if (_mesh_parent->GetBdrAttribute(i) == _elec_attrs.first)
    {
      bdr_els.push_back(i);
    }
  }

  Plane3D plane;

  if (bdr_els.size() > 0)
  {
    plane.Make3DPlane(_mesh_parent, _mesh_parent->GetBdrElementFaceIndex(bdr_els[0]));
  }

  std::vector<int> elec_vtx;
  // Create a vector containing all of the vertices on the electrode
  for (auto b_fc : bdr_els)
  {

    mfem::Array<int> face_vtx;
    _mesh_parent->GetFaceVertices(_mesh_parent->GetBdrElementFaceIndex(b_fc), face_vtx);

    for (auto v : face_vtx)
      pushIfUnique(elec_vtx, v);
  }

  // Now we need to find all elements in the mesh that touch, on at least one
  // vertex, the electrode face if they do touch the vertex, are on one side of
  // the electrode, and belong to the coil domain, we add them to our wedge

  std::vector<int> wedge_els;

  for (int e = 0; e < _mesh_parent->GetNE(); ++e)
  {

    if (!IsInDomain(e, _coil_domains, _mesh_parent) ||
        plane.Side(ElementCentre(e, _mesh_parent)) == 1)
      continue;

    mfem::Array<int> elem_vtx;
    _mesh_parent->GetElementVertices(e, elem_vtx);

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
    _mesh_parent->GetElementFaces(e, el_faces, ori);
    for (auto f : el_faces)
      pushIfUnique(wedge_faces, f);
  }

  for (auto wf : wedge_faces)
  {

    int e1, e2;
    _mesh_parent->GetFaceElements(wf, &e1, &e2);

    // If the face is a coil boundary
    if (!(IsInDomain(e1, _coil_domains, _mesh_parent) &&
          IsInDomain(e2, _coil_domains, _mesh_parent)))
    {
      continue;
    }

    // If the face is not true interior
    if (!(_mesh_parent->FaceIsInterior(wf) ||
          (_mesh_parent->GetFaceInformation(wf).tag == mfem::Mesh::FaceInfoTag::SharedConforming ||
           _mesh_parent->GetFaceInformation(wf).tag ==
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
      if (wf == _mesh_parent->GetBdrElementFaceIndex(b_fc))
      {
        test1 = true;
        break;
      }
    }
    if (test1)
      continue;

    // At last, if the face is none of these things, it must be our second
    // electrode
    auto * new_elem = _mesh_parent->GetFace(wf)->Duplicate(_mesh_parent);
    new_elem->SetAttribute(_elec_attrs.second);
    _mesh_parent->AddBdrElement(new_elem);
  }

  // Only after this do we set the domain attributes
  for (auto e : wedge_els)
    _mesh_parent->SetAttribute(e, _new_domain_attr);

  _transition_domain.Append(_new_domain_attr);
  _coil_domains.Append(_new_domain_attr);

  _mesh_parent->FinalizeTopology();
  _mesh_parent->Finalize();
  _mesh_parent->SetAttributes();
}

void
ClosedCoilSolver::PrepareCoilSubmesh()
{
  _mesh_coil = std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromDomain(*_mesh_parent, _coil_domains));

  _grad_phi_aux_coil_fec =
      std::make_unique<mfem::ND_FECollection>(_order_hcurl, _mesh_coil->Dimension());

  _grad_phi_aux_coil_fes =
      std::make_unique<mfem::ParFiniteElementSpace>(_mesh_coil.get(), _grad_phi_aux_coil_fec.get());

  _grad_phi_aux_coil = std::make_unique<mfem::ParGridFunction>(_grad_phi_aux_coil_fes.get());
  *_grad_phi_aux_coil = 0.0;

  _h1_fe_space_coil_fec =
      std::make_unique<mfem::H1_FECollection>(_order_h1, _mesh_coil->Dimension());

  _h1_fe_space_coil =
      std::make_unique<mfem::ParFiniteElementSpace>(_mesh_coil.get(), _h1_fe_space_coil_fec.get());

  _v_coil = std::make_unique<mfem::ParGridFunction>(_h1_fe_space_coil.get());
  *_v_coil = 0.0;

  _mesh_t = std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromDomain(*_mesh_parent, _transition_domain));
}

void
ClosedCoilSolver::SolveTransition()
{
  auto v_parent = std::make_shared<mfem::ParGridFunction>(_h1_fe_space_parent.get());
  *v_parent = 0.0;

  hephaestus::FESpaces fespaces;
  hephaestus::BCMap bc_maps;

  hephaestus::Coefficients coefs;
  coefs._scalars.Register("electrical_conductivity", _sigma);

  hephaestus::GridFunctions gridfunctions;
  gridfunctions.Register("GradPhi_parent", _grad_phi_parent);
  gridfunctions.Register("V_parent", v_parent);

  hephaestus::InputParameters ocs_params;
  ocs_params.SetParam("GradPotentialName", std::string("GradPhi_parent"));
  ocs_params.SetParam("ConductivityCoefName", std::string("electrical_conductivity"));
  ocs_params.SetParam("IFuncCoefName", std::string("I"));
  ocs_params.SetParam("PotentialName", std::string("V_parent"));
  ocs_params.SetParam("SolverOptions", _solver_options);

  hephaestus::OpenCoilSolver opencoil(ocs_params, _transition_domain, _elec_attrs);

  opencoil.Init(gridfunctions, fespaces, bc_maps, coefs);
  opencoil.Apply(_final_lf.get());

  _mesh_coil->Transfer(*v_parent, *_v_coil);
}

void
ClosedCoilSolver::SolveCoil()
{
  // -(σ∇Va,∇ψ) = (σ∇Vt,∇ψ)
  // where Va is Vaux_coil_, the auxiliary continuous "potential"
  // ψ are the H1 test functions
  // Vt is the transition potential
  // The boundary terms are zero because ∇Va and ∇Vt are perpendicular
  // to the coil boundaries

  mfem::ParGridFunction vaux_coil(_h1_fe_space_coil.get());
  vaux_coil = 0.0;

  mfem::ParBilinearForm a_t(_h1_fe_space_coil.get());
  mfem::ParLinearForm b_coil(_h1_fe_space_coil.get());
  b_coil = 0.0;

  attrToMarker(_transition_domain, _transition_markers, _mesh_coil->attributes.Max());
  a_t.AddDomainIntegrator(new mfem::DiffusionIntegrator(*_sigma), _transition_markers);
  a_t.Assemble();
  a_t.Finalize();
  a_t.AddMult(*_v_coil, b_coil, 1.0);

  mfem::ParBilinearForm a_coil(_h1_fe_space_coil.get());
  a_coil.AddDomainIntegrator(new mfem::DiffusionIntegrator(*_sigma));
  a_coil.Assemble();

  mfem::Array<int> ess_bdr_tdofs_coil;

  // This creates a binary representation of which MPI ranks contain at
  // least one element
  int ref_rank = 0;
  int has_els = (bool)_mesh_coil->GetNE() ? 1 << mfem::Mpi::WorldRank() : 0;
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

  mfem::HypreParMatrix a0_coil;
  mfem::Vector x0_coil;
  mfem::Vector b0_coil;
  a_coil.FormLinearSystem(ess_bdr_tdofs_coil, vaux_coil, b_coil, a0_coil, x0_coil, b0_coil);
  hephaestus::DefaultH1PCGSolver a_coil_solver(_solver_options, a0_coil);
  a_coil_solver.Mult(b0_coil, x0_coil);
  a_coil.RecoverFEMSolution(x0_coil, b_coil, vaux_coil);

  // Now we form the final coil current
  mfem::ParDiscreteLinearOperator grad(_h1_fe_space_coil.get(), _grad_phi_aux_coil->ParFESpace());
  grad.AddDomainInterpolator(new mfem::GradientInterpolator());
  grad.Assemble();
  grad.Mult(vaux_coil, *_grad_phi_aux_coil);

  if (_grad_phi_transfer)
    _grad_phi_t_parent = std::make_unique<mfem::ParGridFunction>(*_grad_phi_parent);

  *_grad_phi_parent = 0.0;
  _mesh_coil->Transfer(*_grad_phi_aux_coil, *_grad_phi_parent);

  mfem::ParBilinearForm m1(_h_curl_fe_space_parent.get());
  hephaestus::attrToMarker(_coil_domains, _coil_markers, _mesh_parent->attributes.Max());

  m1.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(_sigma.get()), _coil_markers);
  m1.Assemble();
  m1.AddMult(*_grad_phi_parent, *_final_lf, -1.0);

  // We can't properly calculate the flux of Jaux on the parent mesh, so we
  // transfer it first to the transition mesh. This will be used in the
  // normalisation step
  auto grad_phi_aux_t_fec =
      std::make_unique<mfem::ND_FECollection>(_order_hcurl, _mesh_t->Dimension());

  auto grad_phi_aux_t_pfes =
      std::make_unique<mfem::ParFiniteElementSpace>(_mesh_t.get(), grad_phi_aux_t_fec.get());

  auto grad_phi_aux_t = std::make_unique<mfem::ParGridFunction>(grad_phi_aux_t_pfes.get());
  *grad_phi_aux_t = 0.0;

  _mesh_t->Transfer(*_grad_phi_parent, *grad_phi_aux_t);

  // The total flux across the electrode face is Φ_t-Φ_aux
  // where Φ_t is the transition flux, already normalised to be 1
  double flux = 1.0 - calcFlux(grad_phi_aux_t.get(), _elec_attrs.first, *_sigma);

  if (_grad_phi_transfer)
  {
    *_grad_phi_parent -= *_grad_phi_t_parent;
    *_grad_phi_parent /= -flux;
    *_grad_phi_t_parent = *_grad_phi_parent;
  }

  *_final_lf /= flux;
}

void
ClosedCoilSolver::RestoreAttributes()
{
  // Domain attributes
  for (int e = 0; e < _mesh_parent->GetNE(); ++e)
  {
    _mesh_parent->SetAttribute(e, _old_dom_attrs[e]);
  }

  _mesh_parent->FinalizeTopology();
  _mesh_parent->Finalize();
  _mesh_parent->SetAttributes();
}

// Auxiliary methods

bool
ClosedCoilSolver::IsInDomain(const int el, const mfem::Array<int> & dom, const mfem::ParMesh * mesh)
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
ClosedCoilSolver::IsInDomain(const int el, const int & sd, const mfem::ParMesh * mesh)
{

  // This is for ghost elements
  if (el < 0)
    return false;

  return mesh->GetAttribute(el) == sd;
}

mfem::Vector
ClosedCoilSolver::ElementCentre(int el, mfem::ParMesh * pm)
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

Plane3D::Plane3D()
{
  _u = std::make_unique<mfem::Vector>(3);
  *_u = 0.0;
}

void
Plane3D::Make3DPlane(const mfem::ParMesh * pm, const int face)
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
  v[0].cross3D(v[1], *_u);
  *_u /= _u->Norml2();

  // Finally, we find d:
  _d = *_u * v[2];
}

int
Plane3D::Side(const mfem::Vector v)
{
  double val = *_u * v - _d;

  if (val > 0)
    return 1;
  else if (val < 0)
    return -1;
  else
    return 0;
}

} // namespace hephaestus