// Solves:
// ∇×(ν∇×A) + σ(dA/dt + ∇ V) = Jᵉ
// ∇·(σ(dA/dt + ∇ V))= 0

// where
// Jᵉ ∈ H(div) source field
// A ∈ H(curl)
// V ∈ H1

//* in weak form
//* (ν∇×A, ∇×A') + (σ(dA/dt + ∇ V), A') - (Jᵉ, A') - <(ν∇×A) × n, A'>  = 0
//* (σ(dA/dt + ∇ V), ∇V') - <σ(dA/dt + ∇ V)·n, V'> =0
//*
//* where:
//* reluctivity ν = 1/μ
//* electrical_conductivity σ=1/ρ
//* Magnetic vector potential A
//* Scalar electric potential V
//* Electric field, E = -dA/dt -∇V
//* Magnetic flux density, B = ∇×A
//* Magnetic field H = ν∇×A
//* Current density J = -σ(dA/dt + ∇ V)

//* Either:
//* B.n (or E×n) at boundary: A×n (Dirichlet)
//* H×n at boundary: ν∇×A (Integrated)
//* -σ(dA/dt + ∇ V)·n (J·n, Neumann), V (potential, Dirichlet)

#include "av_formulation.hpp"

#include <utility>

namespace hephaestus
{

AVFormulation::AVFormulation(std::string alpha_coef_name,
                             std::string inv_alpha_coef_name,
                             std::string beta_coef_name,
                             std::string vector_potential_name,
                             std::string scalar_potential_name)
  : _alpha_coef_name(std::move(alpha_coef_name)),
    _inv_alpha_coef_name(std::move(inv_alpha_coef_name)),
    _beta_coef_name(std::move(beta_coef_name)),
    _vector_potential_name(std::move(vector_potential_name)),
    _scalar_potential_name(std::move(scalar_potential_name))
{
}

void
AVFormulation::ConstructEquationSystem()
{
  hephaestus::InputParameters av_system_params;
  av_system_params.SetParam("VectorPotentialName", _vector_potential_name);
  av_system_params.SetParam("ScalarPotentialName", _scalar_potential_name);
  av_system_params.SetParam("AlphaCoefName", _alpha_coef_name);
  av_system_params.SetParam("BetaCoefName", _beta_coef_name);
  GetProblem()->_td_equation_system =
      std::make_unique<hephaestus::AVEquationSystem>(av_system_params);
}

void
AVFormulation::ConstructOperator()
{
  _problem->_td_operator = std::make_unique<hephaestus::AVOperator>(*(_problem->_pmesh),
                                                                    _problem->_fespaces,
                                                                    _problem->_gridfunctions,
                                                                    _problem->_bc_map,
                                                                    _problem->_coefficients,
                                                                    _problem->_sources,
                                                                    _problem->_solver_options);
  _problem->_td_operator->SetEquationSystem(_problem->_td_equation_system.get());
  _problem->_td_operator->SetGridFunctions();
};

void
AVFormulation::RegisterGridFunctions()
{
  int & myid = GetProblem()->_myid;
  hephaestus::GridFunctions & gridfunctions = GetProblem()->_gridfunctions;

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(_vector_potential_name))
  {
    if (myid == 0)
    {
      MFEM_WARNING(_vector_potential_name
                   << " not found in gridfunctions: building gridfunction from "
                      "defaults");
    }
    AddFESpace(std::string("_HCurlFESpace"), std::string("ND_3D_P2"));
    AddGridFunction(_vector_potential_name, std::string("_HCurlFESpace"));
  }

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(_scalar_potential_name))
  {
    if (myid == 0)
    {
      MFEM_WARNING(_scalar_potential_name
                   << " not found in gridfunctions: building gridfunction from "
                      "defaults");
    }
    AddFESpace(std::string("_H1FESpace"), std::string("H1_3D_P2"));
    AddGridFunction(_scalar_potential_name, std::string("_H1FESpace"));
  }

  // Register time derivatives
  TimeDomainProblemBuilder::RegisterGridFunctions();
};

void
AVFormulation::RegisterCoefficients()
{
  hephaestus::Coefficients & coefficients = GetProblem()->_coefficients;
  if (!coefficients._scalars.Has(_inv_alpha_coef_name))
  {
    MFEM_ABORT(_inv_alpha_coef_name + " coefficient not found.");
  }
  if (!coefficients._scalars.Has(_beta_coef_name))
  {
    MFEM_ABORT(_beta_coef_name + " coefficient not found.");
  }

  coefficients._scalars.Register(
      _alpha_coef_name,
      new mfem::TransformedCoefficient(
          &_one_coef, coefficients._scalars.Get(_inv_alpha_coef_name), fracFunc),
      true);
}

AVEquationSystem::AVEquationSystem(const hephaestus::InputParameters & params)
  : TimeDependentEquationSystem(params),
    _a_name(params.GetParam<std::string>("VectorPotentialName")),
    _v_name(params.GetParam<std::string>("ScalarPotentialName")),
    _alpha_coef_name(params.GetParam<std::string>("AlphaCoefName")),
    _beta_coef_name(params.GetParam<std::string>("BetaCoefName")),
    _dtalpha_coef_name(std::string("dt_") + _alpha_coef_name),
    _neg_beta_coef_name(std::string("negative_") + _beta_coef_name),
    _neg_coef(-1.0)
{
}

void
AVEquationSystem::Init(hephaestus::GridFunctions & gridfunctions,
                       const hephaestus::FESpaces & fespaces,
                       hephaestus::BCMap & bc_map,
                       hephaestus::Coefficients & coefficients)
{
  coefficients._scalars.Register(
      _dtalpha_coef_name,
      new mfem::TransformedCoefficient(
          &_dt_coef, coefficients._scalars.Get(_alpha_coef_name), prodFunc),
      true);

  coefficients._scalars.Register(
      _neg_beta_coef_name,
      new mfem::TransformedCoefficient(
          &_neg_coef, coefficients._scalars.Get(_beta_coef_name), prodFunc),
      true);

  TimeDependentEquationSystem::Init(gridfunctions, fespaces, bc_map, coefficients);
}

void
AVEquationSystem::AddKernels()
{
  AddVariableNameIfMissing(_a_name);
  std::string da_dt_name = GetTimeDerivativeName(_a_name);
  AddVariableNameIfMissing(_v_name);
  std::string dv_dt_name = GetTimeDerivativeName(_v_name);

  // (α∇×A_{n}, ∇×A')
  hephaestus::InputParameters weak_curl_curl_params;
  weak_curl_curl_params.SetParam("CoupledVariableName", _a_name);
  weak_curl_curl_params.SetParam("CoefficientName", _alpha_coef_name);
  AddKernel(da_dt_name, std::make_shared<hephaestus::WeakCurlCurlKernel>(weak_curl_curl_params));

  // (αdt∇×dA/dt_{n+1}, ∇×A')
  hephaestus::InputParameters curl_curl_params;
  curl_curl_params.SetParam("CoefficientName", _dtalpha_coef_name);
  AddKernel(da_dt_name, std::make_shared<hephaestus::CurlCurlKernel>(curl_curl_params));

  // (βdA/dt_{n+1}, A')
  hephaestus::InputParameters vector_fe_mass_params;
  vector_fe_mass_params.SetParam("CoefficientName", _beta_coef_name);
  AddKernel(da_dt_name, std::make_shared<hephaestus::VectorFEMassKernel>(vector_fe_mass_params));

  // (σ ∇ V, dA'/dt)
  hephaestus::InputParameters mixed_vector_gradient_params;
  mixed_vector_gradient_params.SetParam("CoefficientName", _beta_coef_name);
  AddKernel(_v_name,
            da_dt_name,
            std::make_shared<hephaestus::MixedVectorGradientKernel>(mixed_vector_gradient_params));

  // (σ ∇ V, ∇ V')
  hephaestus::InputParameters diffusion_params;
  diffusion_params.SetParam("CoefficientName", _beta_coef_name);
  AddKernel(_v_name, std::make_shared<hephaestus::DiffusionKernel>(diffusion_params));

  // (σdA/dt, ∇ V')
  hephaestus::InputParameters vector_fe_weak_divergence_params;
  vector_fe_weak_divergence_params.SetParam("CoefficientName", _beta_coef_name);
  AddKernel(
      da_dt_name,
      _v_name,
      std::make_shared<hephaestus::VectorFEWeakDivergenceKernel>(vector_fe_weak_divergence_params));
}

AVOperator::AVOperator(mfem::ParMesh & pmesh,
                       hephaestus::FESpaces & fespaces,
                       hephaestus::GridFunctions & gridfunctions,
                       hephaestus::BCMap & bc_map,
                       hephaestus::Coefficients & coefficients,
                       hephaestus::Sources & sources,
                       hephaestus::InputParameters & solver_options)
  : TimeDomainEquationSystemOperator(
        pmesh, fespaces, gridfunctions, bc_map, coefficients, sources, solver_options)
{
  // Initialize MPI gridfunctions
  MPI_Comm_size(pmesh.GetComm(), &_num_procs);
  MPI_Comm_rank(pmesh.GetComm(), &_myid);

  _aux_var_names.resize(2);
  _aux_var_names.at(0) = "electric_field";
  _aux_var_names.at(1) = "magnetic_flux_density";
}

/*
This is the main computational code that computes dX/dt implicitly
where X is the state vector containing p, u and v.

Unknowns
s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
du/dt_{n+1} ∈ H(curl)
p_{n+1} ∈ H1

Fully discretised equations
-(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
(α∇×u_{n}, ∇×u') + (αdt∇×du/dt_{n+1}, ∇×u') + (βdu/dt_{n+1}, u')
- (s0_{n+1}, u') - <(α∇×u_{n+1}) × n, u'> = 0
using
u_{n+1} = u_{n} + dt du/dt_{n+1}
*/
void
AVOperator::ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt)
{
  for (unsigned int ind = 0; ind < _local_test_vars.size(); ++ind)
  {
    _local_test_vars.at(ind)->MakeRef(
        _local_test_vars.at(ind)->ParFESpace(), const_cast<mfem::Vector &>(X), _true_offsets[ind]);
    _local_trial_vars.at(ind)->MakeRef(
        _local_trial_vars.at(ind)->ParFESpace(), dX_dt, _true_offsets[ind]);
  }
  _coefficients.SetTime(GetTime());
  _equation_system->SetTimeStep(dt);
  _equation_system->UpdateEquationSystem(_bc_map, _sources);

  _equation_system->FormLinearSystem(_block_a, _true_x, _true_rhs);

  _solver = std::make_unique<hephaestus::DefaultGMRESSolver>(_solver_options,
                                                             *_block_a.As<mfem::HypreParMatrix>());
  // solver = new hephaestus::DefaultGMRESSolver(_solver_options, *blockA,
  //                                             pmesh_->GetComm());

  _solver->Mult(_true_rhs, _true_x);
  _equation_system->RecoverFEMSolution(_true_x, _gridfunctions);
}
} // namespace hephaestus
