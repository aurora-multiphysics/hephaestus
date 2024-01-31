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

OpenCoilSolver::OpenCoilSolver(std::string source_efield_gf_name,
                               std::string phi_gf_name,
                               std::string i_coef_name,
                               std::string cond_coef_name,
                               mfem::Array<int> coil_dom,
                               const std::pair<int, int> electrodes,
                               bool grad_phi_transfer,
                               std::string source_jfield_gf_name,
                               hephaestus::InputParameters solver_options)
  : _source_efield_gf_name(std::move(source_efield_gf_name)),
    _source_jfield_gf_name(std::move(source_jfield_gf_name)),
    _phi_gf_name(std::move(phi_gf_name)),
    _i_coef_name(std::move(i_coef_name)),
    _cond_coef_name(std::move(cond_coef_name)),
    _grad_phi_transfer(std::move(grad_phi_transfer)),
    _solver_options(std::move(solver_options)),
    _coil_domains(std::move(coil_dom)),
    _elec_attrs(electrodes),
    _high_src(highV),
    _low_src(lowV),
    _high_terminal(1),
    _low_terminal(1),
    _ref_face(_elec_attrs.first)
{
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
    std::cout << "Warning: source electric field undefined. The GridFunction "
                 "associated with it will be set to zero.\n";

    _sigma = new mfem::ConstantCoefficient(1.0);
    _owns_sigma = true;

    _grad_phi_transfer = false;
  }

  _source_electric_field = gridfunctions.Get(_source_efield_gf_name);
  if (_source_electric_field == nullptr)
  {
    const std::string error_message = _source_efield_gf_name + " not found in gridfunctions when "
                                                               "creating OpenCoilSolver\n";
    mfem::mfem_error(error_message.c_str());
  }
  else if (_source_electric_field->ParFESpace()->FEColl()->GetContType() !=
           mfem::FiniteElementCollection::TANGENTIAL)
  {
    mfem::mfem_error("Electric field GridFunction must be of HCurl type.");
  }

  if (!_source_jfield_gf_name.empty())
  {
    _source_current_density = gridfunctions.Get(_source_jfield_gf_name);
    if (_source_current_density == nullptr)
    {
      const std::string error_message = _source_jfield_gf_name + " not found in gridfunctions when "
                                                                 "creating OpenCoilSolver\n";
      mfem::mfem_error(error_message.c_str());
    }
    else if (_source_current_density->ParFESpace()->FEColl()->GetContType() !=
             mfem::FiniteElementCollection::NORMAL)
    {
      mfem::mfem_error("Current density GridFunction must be of HDiv type.");
    }
  }

  _order_hcurl = _source_electric_field->ParFESpace()->FEColl()->GetOrder();

  _phi_parent = gridfunctions.Get(_phi_gf_name);
  if (_phi_parent == nullptr)
  {
    std::cout << _phi_gf_name + " not found in gridfunctions when "
                                "creating OpenCoilSolver.\n";
    _order_h1 = _order_hcurl;
  }
  else if (_phi_parent->ParFESpace()->FEColl()->GetContType() !=
           mfem::FiniteElementCollection::CONTINUOUS)
  {
    mfem::mfem_error("V GridFunction must be of H1 type.");
  }
  else
  {
    _order_h1 = _phi_parent->ParFESpace()->FEColl()->GetOrder();
    _phi_t_parent = std::make_unique<mfem::ParGridFunction>(*_phi_parent);
  }

  _mesh_parent = _source_electric_field->ParFESpace()->GetParMesh();

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
      mfem::IntRules.Get(_source_electric_field->ParFESpace()->GetFE(0)->GetGeomType(), 1)
          .IntPoint(0);

  double i = _itotal->Eval(*tr, ip);

  if (_grad_phi_transfer)
  {
    *_source_electric_field = 0.0;
    _source_electric_field->Add(i, *_grad_phi_t_parent);
  }

  if (_phi_parent != nullptr)
  {
    *_phi_parent = 0.0;
    _phi_parent->Add(i, *_phi_t_parent);
  }

  if (_source_current_density)
  {
    hephaestus::GridFunctions aux_gf;
    aux_gf.Register("source_electric_field", _source_electric_field, false);
    aux_gf.Register("source_current_density", _source_current_density, false);

    hephaestus::Coefficients aux_coef;
    aux_coef._scalars.Register("electrical_conductivity", _sigma, false);

    hephaestus::ScaledVectorGridFunctionAux current_density_auxsolver(
        "source_electric_field", "source_current_density", "electrical_conductivity", -1.0);
    current_density_auxsolver.Init(aux_gf, aux_coef);
    current_density_auxsolver.Solve();
  }

  lf->Add(i, *_final_lf);
}

void
OpenCoilSolver::InitChildMesh()
{
  if (_mesh_child == nullptr)
    _mesh_child = std::make_unique<mfem::ParSubMesh>(
        mfem::ParSubMesh::CreateFromDomain(*_mesh_parent, _coil_domains));
}

void
OpenCoilSolver::MakeFESpaces()
{
  if (_h1_fe_space_child == nullptr)
  {
    _h1_fe_space_fec_child =
        std::make_unique<mfem::H1_FECollection>(_order_h1, _mesh_child->Dimension());
    _h1_fe_space_child = std::make_unique<mfem::ParFiniteElementSpace>(
        _mesh_child.get(), _h1_fe_space_fec_child.get());
  }

  if (_h_curl_fe_space_child == nullptr)
  {
    _h_curl_fe_space_fec_child =
        std::make_unique<mfem::ND_FECollection>(_order_hcurl, _mesh_child->Dimension());
    _h_curl_fe_space_child = std::make_unique<mfem::ParFiniteElementSpace>(
        _mesh_child.get(), _h_curl_fe_space_fec_child.get());
  }
}

void
OpenCoilSolver::MakeGridFunctions()
{

  if (_phi_child == nullptr)
    _phi_child = std::make_unique<mfem::ParGridFunction>(_h1_fe_space_child.get());

  if (_grad_phi_child == nullptr)
    _grad_phi_child = std::make_unique<mfem::ParGridFunction>(_h_curl_fe_space_child.get());

  if (_grad_phi_t_parent == nullptr)
    _grad_phi_t_parent = std::make_unique<mfem::ParGridFunction>(*_source_electric_field);

  *_phi_child = 0.0;
  *_grad_phi_child = 0.0;
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
                    new hephaestus::ScalarDirichletBC(std::string("V"), _high_terminal, &_high_src),
                    true);

  _bc_maps.Register("low_potential",
                    new hephaestus::ScalarDirichletBC(std::string("V"), _low_terminal, &_low_src),
                    true);

  // NB: register false to avoid double-free.
  hephaestus::FESpaces fespaces;
  fespaces.Register(std::string("HCurl"), _h_curl_fe_space_child.get(), false);
  fespaces.Register(std::string("H1"), _h1_fe_space_child.get(), false);

  // NB: register false to avoid double-free.
  hephaestus::GridFunctions gridfunctions;
  gridfunctions.Register(std::string("GradPhi"), _grad_phi_child.get(), false);
  gridfunctions.Register(std::string("V"), _phi_child.get(), false);

  hephaestus::Coefficients coefs;
  coefs._scalars.Register("electric_conductivity", _sigma, false);

  hephaestus::ScalarPotentialSource sps(
      "GradPhi", "V", "HCurl", "H1", "electric_conductivity", _solver_options);
  sps.Init(gridfunctions, fespaces, _bc_maps, coefs);

  mfem::ParLinearForm dummy(_h_curl_fe_space_child.get());
  sps.Apply(&dummy);

  // Normalise the current through the wedges and use them as a reference
  double flux = calcFlux(_grad_phi_child.get(), _ref_face, *_sigma);
  *_grad_phi_child /= abs(flux);
  if (_phi_child)
    *_phi_child /= abs(flux);

  _mesh_child->Transfer(*_grad_phi_child, *_grad_phi_t_parent);
  if (_phi_parent)
    _mesh_child->Transfer(*_phi_child, *_phi_t_parent);

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
    _m1 = std::make_unique<mfem::ParBilinearForm>(_source_electric_field->ParFESpace());
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