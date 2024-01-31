#include "scalar_potential_source.hpp"

namespace hephaestus
{

// Create a scalar potential source that can add terms of the form
// J0_{n+1} ∈ H(div) source field, where J0 = -β∇p and β is a conductivity
// coefficient.
ScalarPotentialSource::ScalarPotentialSource(std::string grad_phi_gf_name,
                                             std::string phi_gf_name,
                                             std::string hcurl_fespace_name,
                                             std::string h1_fespace_name,
                                             std::string coef_name,
                                             hephaestus::InputParameters solver_options)
  : _grad_phi_gf_name(std::move(grad_phi_gf_name)),
    _phi_gf_name(std::move(phi_gf_name)),
    _hcurl_fespace_name(std::move(hcurl_fespace_name)),
    _h1_fespace_name(std::move(h1_fespace_name)),
    _coef_name(std::move(coef_name)),
    _solver_options(std::move(solver_options))
{
}

void
ScalarPotentialSource::Init(hephaestus::GridFunctions & gridfunctions,
                            const hephaestus::FESpaces & fespaces,
                            hephaestus::BCMap & bc_map,
                            hephaestus::Coefficients & coefficients)
{
  _h1_fe_space = fespaces.Get(_h1_fespace_name);
  if (_h1_fe_space == nullptr)
  {
    const std::string error_message = _h1_fespace_name + " not found in fespaces when "
                                                         "creating ScalarPotentialSource\n";
    mfem::mfem_error(error_message.c_str());
  }

  _h_curl_fe_space = fespaces.Get(_hcurl_fespace_name);
  if (_h_curl_fe_space == nullptr)
  {
    const std::string error_message = _hcurl_fespace_name + " not found in fespaces when "
                                                            "creating ScalarPotentialSource\n";
    mfem::mfem_error(error_message.c_str());
  }

  _phi = gridfunctions.Get(_phi_gf_name);
  if (_phi == nullptr)
  {
    std::cout << _phi_gf_name + " not found in gridfunctions when "
                                "creating ScalarPotentialSource. "
                                "Creating new ParGridFunction\n";
    _phi = new mfem::ParGridFunction(_h1_fe_space);
    gridfunctions.Register(_phi_gf_name, _phi, true);
  }

  _grad_phi = gridfunctions.Get(_grad_phi_gf_name);
  if (_grad_phi == nullptr)
  {
    _grad_phi = new mfem::ParGridFunction(_h_curl_fe_space);
    gridfunctions.Register(_grad_phi_gf_name, _grad_phi, true);
  }

  _bc_map = &bc_map;

  _beta_coef = coefficients._scalars.Get(_coef_name);

  _a0 = std::make_unique<mfem::ParBilinearForm>(_h1_fe_space);
  _a0->AddDomainIntegrator(new mfem::DiffusionIntegrator(*_beta_coef));
  _a0->Assemble();

  BuildGrad();
  BuildM1(_beta_coef);
  // a0(p, p') = (β ∇ p, ∇ p')

  _b0 = std::make_unique<mfem::ParLinearForm>(_h1_fe_space);
  _diffusion_mat = std::make_unique<mfem::HypreParMatrix>();
  _p_tdofs = std::make_unique<mfem::Vector>();
  _b0_tdofs = std::make_unique<mfem::Vector>();
}

ScalarPotentialSource::~ScalarPotentialSource() = default;

void
ScalarPotentialSource::BuildM1(mfem::Coefficient * Sigma)
{
  _m1 = std::make_unique<mfem::ParBilinearForm>(_h_curl_fe_space);
  _m1->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*Sigma));
  _m1->Assemble();

  // Don't finalize or parallel assemble this is done in FormLinearSystem.
}

void
ScalarPotentialSource::BuildGrad()
{
  _grad = std::make_unique<mfem::ParDiscreteLinearOperator>(_h1_fe_space, _h_curl_fe_space);
  _grad->AddDomainInterpolator(new mfem::GradientInterpolator());
  _grad->Assemble();

  // no ParallelAssemble since this will be applied to GridFunctions
}

void
ScalarPotentialSource::Apply(mfem::ParLinearForm * lf)
{
  // -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
  // a0(p_{n+1}, p') = b0(p')
  // a0(p, p') = (β ∇ p, ∇ p')
  // b0(p') = <n.s0, p'>
  mfem::ParGridFunction phi_gf(_h1_fe_space);
  mfem::Array<int> poisson_ess_tdof_list;
  phi_gf = 0.0;
  *_b0 = 0.0;
  _bc_map->ApplyEssentialBCs(
      _phi_gf_name, poisson_ess_tdof_list, phi_gf, (_h1_fe_space->GetParMesh()));
  _bc_map->ApplyIntegratedBCs(_phi_gf_name, *_b0, (_h1_fe_space->GetParMesh()));
  _b0->Assemble();

  _a0->Update();
  _a0->Assemble();
  _a0->FormLinearSystem(
      poisson_ess_tdof_list, phi_gf, *_b0, *_diffusion_mat, *_p_tdofs, *_b0_tdofs);

  if (_a0_solver == nullptr)
  {
    _a0_solver = std::make_unique<hephaestus::DefaultH1PCGSolver>(_solver_options, *_diffusion_mat);
  }
  // Solve
  _a0_solver->Mult(*_b0_tdofs, *_p_tdofs);

  // "undo" the static condensation saving result in grid function dP
  _a0->RecoverFEMSolution(*_p_tdofs, *_b0, *_phi);

  _grad->Mult(*_phi, *_grad_phi);

  _m1->Update();
  _m1->Assemble();
  _m1->AddMult(*_grad_phi, *lf, 1.0);
}

void
ScalarPotentialSource::SubtractSource(mfem::ParGridFunction * gf)
{
  _grad->AddMult(*_phi, *gf, -1.0);
}

} // namespace hephaestus
