#include "open_coil.hpp"

#include <utility>

namespace hephaestus
{

///// THESE FUNCTIONS WILL EVENTUALLY GO INTO A UTILS FILE ///////////

double
highV(const mfem::Vector & x, double t)
{
  return 0.5;
}
double
lowV(const mfem::Vector & x, double t)
{
  return -0.5;
}

double
calcFlux(mfem::GridFunction * v_field, int face_attr, mfem::Coefficient & q)
{

  double flux = 0.0;
  double area = 0.0;

  mfem::FiniteElementSpace * fes = v_field->FESpace();
  mfem::Mesh * mesh = fes->GetMesh();

  mfem::Vector local_dofs, normal_vec;
  mfem::DenseMatrix dshape;
  mfem::Array<int> dof_ids;

  for (int i = 0; i < mesh->GetNBE(); i++)
  {

    if (mesh->GetBdrAttribute(i) != face_attr)
      continue;

    mfem::FaceElementTransformations * f_tr =
        mesh->GetFaceElementTransformations(mesh->GetBdrElementFaceIndex(i));
    if (f_tr == nullptr)
      continue;

    const mfem::FiniteElement & elem = *fes->GetFE(f_tr->Elem1No);
    const int int_order = 2 * elem.GetOrder() + 3;
    const mfem::IntegrationRule & ir = mfem::IntRules.Get(f_tr->FaceGeom, int_order);

    fes->GetElementDofs(f_tr->Elem1No, dof_ids);
    v_field->GetSubVector(dof_ids, local_dofs);
    const int space_dim = f_tr->Face->GetSpaceDim();
    normal_vec.SetSize(space_dim);
    dshape.SetSize(elem.GetDof(), space_dim);

    for (int j = 0; j < ir.GetNPoints(); j++)
    {

      const mfem::IntegrationPoint & ip = ir.IntPoint(j);
      mfem::IntegrationPoint eip;
      f_tr->Loc1.Transform(ip, eip);
      f_tr->Face->SetIntPoint(&ip);
      double face_weight = f_tr->Face->Weight();
      double val = 0.0;
      f_tr->Elem1->SetIntPoint(&eip);
      elem.CalcVShape(*f_tr->Elem1, dshape);
      mfem::CalcOrtho(f_tr->Face->Jacobian(), normal_vec);
      val += dshape.InnerProduct(normal_vec, local_dofs) / face_weight;

      // Measure the area of the boundary
      area += ip.weight * face_weight;

      // Integrate alpha * n.Grad(x) + beta * x
      flux += q.Eval(*f_tr, ip) * val * ip.weight * face_weight;
    }
  }

  double total_flux;
  MPI_Allreduce(&flux, &total_flux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return total_flux;
}

double
calcFlux(mfem::GridFunction * v_field, int face_attr)
{
  mfem::ConstantCoefficient one_coef(1.0);
  return calcFlux(v_field, face_attr, one_coef);
}

void
SubdomainToArray(const std::vector<hephaestus::Subdomain> & sd, mfem::Array<int> & arr)
{
  arr.DeleteAll();
  for (auto s : sd)
    arr.Append(s._id);
}

void
SubdomainToArray(const hephaestus::Subdomain & sd, mfem::Array<int> & arr)
{
  arr.DeleteAll();
  arr.Append(sd._id);
}

void
inheritBdrAttributes(const mfem::ParMesh * parent_mesh, mfem::ParSubMesh * child_mesh)
{

  int face, ori, att;
  auto map = child_mesh->GetParentToSubMeshFaceIDMap();

  for (int bdr = 0; bdr < parent_mesh->GetNBE(); ++bdr)
  {

    parent_mesh->GetBdrElementFace(bdr, &face, &ori);
    if (map[face] != -1)
    {
      att = parent_mesh->GetBdrAttribute(bdr);
      auto * new_elem = child_mesh->GetFace(map[face])->Duplicate(child_mesh);
      new_elem->SetAttribute(att);
      child_mesh->AddBdrElement(new_elem);
    }
  }

  child_mesh->FinalizeTopology();
  child_mesh->Finalize();
  child_mesh->SetAttributes();
}

void
attrToMarker(const mfem::Array<int> attr_list, mfem::Array<int> & marker_list, int max_attr)
{

  marker_list.SetSize(max_attr);
  marker_list = 0;

  for (auto a : attr_list)
    marker_list[a - 1] = 1;
}

void
cleanDivergence(mfem::ParGridFunction & Vec_GF, hephaestus::InputParameters solve_pars)
{

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

void
cleanDivergence(hephaestus::GridFunctions & gfs,
                hephaestus::BCMap & bcs,
                const std::string vec_gf_name,
                const std::string scalar_gf_name,
                hephaestus::InputParameters solve_pars)
{

  hephaestus::InputParameters pars;
  hephaestus::FESpaces fes;

  pars.SetParam("VectorGridFunctionName", vec_gf_name);
  pars.SetParam("ScalarGridFunctionName", scalar_gf_name);
  pars.SetParam("SolverOptions", solve_pars);
  hephaestus::HelmholtzProjector projector(pars);
  projector.Project(gfs, fes, bcs);
}

/////////////////////////////////////////////////////////////////////

OpenCoilSolver::OpenCoilSolver(const hephaestus::InputParameters & params,
                               mfem::Array<int> coil_dom,
                               const std::pair<int, int> electrodes)
  : _grad_phi_name(params.GetParam<std::string>("GradPotentialName")),
    _v_gf_name(params.GetParam<std::string>("PotentialName")),
    _i_coef_name(params.GetParam<std::string>("IFuncCoefName")),
    _cond_coef_name(params.GetParam<std::string>("ConductivityCoefName")),
    _grad_phi_transfer(params.GetOptionalParam<bool>("GradPhiTransfer", true)),
    _coil_domains(std::move(coil_dom)),
    _elec_attrs(electrodes),
    _high_src(std::make_shared<mfem::FunctionCoefficient>(highV)),
    _low_src(std::make_shared<mfem::FunctionCoefficient>(lowV)),
    _high_terminal(1),
    _low_terminal(1)
{

  hephaestus::InputParameters default_pars;
  default_pars.SetParam("Tolerance", float(1.0e-20));
  default_pars.SetParam("AbsTolerance", float(1.0e-20));
  default_pars.SetParam("MaxIter", (unsigned int)1000);
  default_pars.SetParam("PrintLevel", 1);

  _solver_options =
      params.GetOptionalParam<hephaestus::InputParameters>("SolverOptions", default_pars);

  _ref_face = _elec_attrs.first;
}

OpenCoilSolver::~OpenCoilSolver()
{
  if (_owns_sigma)
    delete _sigma;
  if (_owns_itotal)
    delete _itotal;
}

void
OpenCoilSolver::Init(hephaestus::GridFunctions & gridfunctions,
                     const hephaestus::FESpaces & fespaces,
                     hephaestus::BCMap & bc_map,
                     hephaestus::Coefficients & coefficients)
{

  _itotal = coefficients._scalars.Get(_i_coef_name);
  if (_itotal == nullptr)
  {
    std::cout << _i_coef_name + " not found in coefficients when "
                                "creating OpenCoilSolver. "
                                "Assuming unit current.\n";
    _itotal = new mfem::ConstantCoefficient(1.0);
    _owns_itotal = true;
  }

  _sigma = coefficients._scalars.Get(_cond_coef_name);
  if (_sigma == nullptr)
  {
    std::cout << _cond_coef_name + " not found in coefficients when "
                                   "creating OpenCoilSolver. "
                                   "Assuming unit conductivity.\n";
    std::cout << "Warning: GradPhi field undefined. The GridFunction "
                 "associated with it will be set to zero.\n";

    _sigma = new mfem::ConstantCoefficient(1.0);
    _owns_sigma = true;

    _grad_phi_transfer = false;
  }

  _grad_phi_parent = gridfunctions.Get(_grad_phi_name);
  if (_grad_phi_parent == nullptr)
  {
    const std::string error_message = _grad_phi_name + " not found in gridfunctions when "
                                                       "creating OpenCoilSolver\n";
    mfem::mfem_error(error_message.c_str());
  }
  else if (_grad_phi_parent->ParFESpace()->FEColl()->GetContType() !=
           mfem::FiniteElementCollection::TANGENTIAL)
  {
    mfem::mfem_error("GradPhi GridFunction must be of HCurl type.");
  }

  _order_hcurl = _grad_phi_parent->ParFESpace()->FEColl()->GetOrder();

  _v_parent = gridfunctions.Get(_v_gf_name);
  if (_v_parent == nullptr)
  {
    std::cout << _v_gf_name + " not found in gridfunctions when "
                              "creating OpenCoilSolver.\n";
    _order_h1 = _order_hcurl;
  }
  else if (_v_parent->ParFESpace()->FEColl()->GetContType() !=
           mfem::FiniteElementCollection::CONTINUOUS)
  {
    mfem::mfem_error("V GridFunction must be of H1 type.");
  }
  else
  {
    _order_h1 = _v_parent->ParFESpace()->FEColl()->GetOrder();
    _vt_parent = std::make_unique<mfem::ParGridFunction>(*_v_parent);
  }

  _mesh_parent = _grad_phi_parent->ParFESpace()->GetParMesh();

  InitChildMesh();
  MakeFESpaces();
  MakeGridFunctions();
  SetBCs();
  SPSCurrent();
}

void
OpenCoilSolver::Apply(mfem::ParLinearForm * lf)
{

  // The transformation and integration points themselves are not relevant, it's
  // just so we can call Eval
  mfem::ElementTransformation * tr = _mesh_parent->GetElementTransformation(0);
  const mfem::IntegrationPoint & ip =
      mfem::IntRules.Get(_grad_phi_parent->ParFESpace()->GetFE(0)->GetGeomType(), 1).IntPoint(0);

  double i = _itotal->Eval(*tr, ip);

  *_grad_phi_parent = 0.0;
  if (_grad_phi_transfer)
    _grad_phi_parent->Add(i, *_grad_phi_t_parent);

  if (_v_parent != nullptr)
  {
    *_v_parent = 0.0;
    _v_parent->Add(i, *_vt_parent);
  }

  lf->Add(i, *_final_lf);
}

void
OpenCoilSolver::SubtractSource(mfem::ParGridFunction * gf)
{
}

void
OpenCoilSolver::InitChildMesh()
{
  if (_mesh == nullptr)
    _mesh = std::make_unique<mfem::ParSubMesh>(
        mfem::ParSubMesh::CreateFromDomain(*_mesh_parent, _coil_domains));
}

void
OpenCoilSolver::MakeFESpaces()
{
  if (_h1_fe_space == nullptr)
  {
    _h1_fe_space_fec = std::make_unique<mfem::H1_FECollection>(_order_h1, _mesh->Dimension());
    _h1_fe_space =
        std::make_unique<mfem::ParFiniteElementSpace>(_mesh.get(), _h1_fe_space_fec.get());
  }

  if (_h_curl_fe_space == nullptr)
  {
    _h_curl_fe_space_fec =
        std::make_unique<mfem::ND_FECollection>(_order_hcurl, _mesh->Dimension());
    _h_curl_fe_space =
        std::make_unique<mfem::ParFiniteElementSpace>(_mesh.get(), _h_curl_fe_space_fec.get());
  }
}

void
OpenCoilSolver::MakeGridFunctions()
{

  if (_v == nullptr)
    _v = std::make_unique<mfem::ParGridFunction>(_h1_fe_space.get());

  if (_grad_phi == nullptr)
    _grad_phi = std::make_unique<mfem::ParGridFunction>(_h_curl_fe_space.get());

  if (_grad_phi_t_parent == nullptr)
    _grad_phi_t_parent = std::make_unique<mfem::ParGridFunction>(*_grad_phi_parent);

  *_v = 0.0;
  *_grad_phi = 0.0;
  *_grad_phi_t_parent = 0.0;
}

void
OpenCoilSolver::SetBCs()
{

  _high_terminal[0] = _elec_attrs.first;
  _low_terminal[0] = _elec_attrs.second;
}

void
OpenCoilSolver::SPSCurrent()
{
  _bc_maps.Register("high_potential",
                    new hephaestus::ScalarDirichletBC(std::string("V"), _high_terminal, _high_src),
                    true);

  _bc_maps.Register("low_potential",
                    new hephaestus::ScalarDirichletBC(std::string("V"), _low_terminal, _low_src),
                    true);

  // NB: register false to avoid double-free.
  hephaestus::FESpaces fespaces;
  fespaces.Register(std::string("HCurl"), _h_curl_fe_space.get(), false);
  fespaces.Register(std::string("H1"), _h1_fe_space.get(), false);

  // NB: register false to avoid double-free.
  hephaestus::GridFunctions gridfunctions;
  gridfunctions.Register(std::string("GradPhi"), _grad_phi.get(), false);
  gridfunctions.Register(std::string("V"), _v.get(), false);

  hephaestus::InputParameters sps_params;
  sps_params.SetParam("GradPotentialName", std::string("GradPhi"));
  sps_params.SetParam("PotentialName", std::string("V"));
  sps_params.SetParam("HCurlFESpaceName", std::string("HCurl"));
  sps_params.SetParam("H1FESpaceName", std::string("H1"));
  sps_params.SetParam("SolverOptions", _solver_options);
  sps_params.SetParam("ConductivityCoefName", std::string("electric_conductivity"));

  hephaestus::Coefficients coefs;
  coefs._scalars.Register("electric_conductivity", _sigma, false);

  hephaestus::ScalarPotentialSource sps(sps_params);
  sps.Init(gridfunctions, fespaces, _bc_maps, coefs);

  mfem::ParLinearForm dummy(_h_curl_fe_space.get());
  sps.Apply(&dummy);

  // Normalise the current through the wedges and use them as a reference
  double flux = calcFlux(_grad_phi.get(), _ref_face, *_sigma);
  *_grad_phi /= abs(flux);
  if (_v)
    *_v /= abs(flux);

  _mesh->Transfer(*_grad_phi, *_grad_phi_t_parent);
  if (_v_parent)
    _mesh->Transfer(*_v, *_vt_parent);

  BuildM1();

  _final_lf = std::make_unique<mfem::ParLinearForm>(_grad_phi_t_parent->ParFESpace());
  *_final_lf = 0.0;
  _m1->AddMult(*_grad_phi_t_parent, *_final_lf, 1.0);
}

void
OpenCoilSolver::BuildM1()
{
  if (_m1 == nullptr)
  {
    _m1 = std::make_unique<mfem::ParBilinearForm>(_grad_phi_parent->ParFESpace());
    hephaestus::attrToMarker(_coil_domains, _coil_markers, _mesh_parent->attributes.Max());
    _m1->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(_sigma), _coil_markers);
    _m1->Assemble();
    _m1->Finalize();
  }
}

void
OpenCoilSolver::SetRefFace(const int face)
{
  _ref_face = face;
}

} // namespace hephaestus