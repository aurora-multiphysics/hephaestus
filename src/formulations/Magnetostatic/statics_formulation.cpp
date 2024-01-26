// Solves the equations
// ∇⋅s0 = 0
// ∇×(α∇×u) = s0

// where
// s0 ∈ H(div) source field
// u ∈ H(curl)

// Dirichlet boundaries constrain u
// Integrated boundaries constrain (α∇×u) × n

// Weak form (Space discretisation)
// (α∇×u, ∇×u') - (s0, u') - <(α∇×u) × n, u'> = 0

// Time discretisation using implicit scheme:
// Unknowns
// u ∈ H(curl)

// Fully discretised equations
// (α∇×u, ∇×u') - (s0, u') - <(α∇×u) × n, u'> = 0

// Rewritten as
// a1(u, u') = b1(u')

// where
// a1(u, u') = (α∇×u, ∇×u')
// b1(u') = (s0, u') + <(α∇×u) × n, u'>
#include "statics_formulation.hpp"

#include <utility>

namespace hephaestus
{

StaticsFormulation::StaticsFormulation(std::string alpha_coef_name, std::string h_curl_var_name)
  : _alpha_coef_name(std::move(alpha_coef_name)), _h_curl_var_name(std::move(h_curl_var_name))
{
}

void
StaticsFormulation::ConstructOperator()
{
  hephaestus::InputParameters & solver_options = GetProblem()->_solver_options;
  solver_options.SetParam("HCurlVarName", _h_curl_var_name);
  solver_options.SetParam("StiffnessCoefName", _alpha_coef_name);
  _problem->_eq_sys_operator =
      std::make_unique<hephaestus::StaticsOperator>(*(_problem->_pmesh),
                                                    _problem->_fespaces,
                                                    _problem->_gridfunctions,
                                                    _problem->_bc_map,
                                                    _problem->_coefficients,
                                                    _problem->_sources,
                                                    _problem->_solver_options);
  _problem->GetOperator()->SetGridFunctions();
};

void
StaticsFormulation::RegisterGridFunctions()
{
  int & myid = GetProblem()->_myid;
  hephaestus::GridFunctions & gridfunctions = GetProblem()->_gridfunctions;
  hephaestus::FESpaces & fespaces = GetProblem()->_fespaces;

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(_h_curl_var_name))
  {
    if (myid == 0)
    {
      MFEM_WARNING(_h_curl_var_name << " not found in gridfunctions: building "
                                       "gridfunction from defaults");
    }
    AddFESpace(std::string("_HCurlFESpace"), std::string("ND_3D_P2"));
    AddGridFunction(_h_curl_var_name, std::string("_HCurlFESpace"));
  };
};

void
StaticsFormulation::RegisterCoefficients()
{
  hephaestus::Coefficients & coefficients = GetProblem()->_coefficients;
  if (!coefficients._scalars.Has(_alpha_coef_name))
  {
    MFEM_ABORT(_alpha_coef_name + " coefficient not found.");
  }
}

StaticsOperator::StaticsOperator(mfem::ParMesh & pmesh,
                                 hephaestus::FESpaces & fespaces,
                                 hephaestus::GridFunctions & gridfunctions,
                                 hephaestus::BCMap & bc_map,
                                 hephaestus::Coefficients & coefficients,
                                 hephaestus::Sources & sources,
                                 hephaestus::InputParameters & solver_options)
  : EquationSystemOperator(
        pmesh, fespaces, gridfunctions, bc_map, coefficients, sources, solver_options),
    _h_curl_var_name(solver_options.GetParam<std::string>("HCurlVarName")),
    _stiffness_coef_name(solver_options.GetParam<std::string>("StiffnessCoefName"))
{
}

void
StaticsOperator::SetGridFunctions()
{
  _state_var_names.push_back(_h_curl_var_name);
  EquationSystemOperator::SetGridFunctions();
};

void
StaticsOperator::Init(mfem::Vector & X)
{
  EquationSystemOperator::Init(X);
  _stiff_coef = _coefficients._scalars.Get(_stiffness_coef_name);
}

/*
This is the main method that solves for u.

Unknowns
s0 ∈ H(div) divergence-free source field
u ∈ H(curl)

Fully discretised equations
(α∇×u, ∇×u') - (s0, u') - <(α∇×u) × n, u'> = 0
*/
void
StaticsOperator::Solve(mfem::Vector & X)
{
  mfem::ParGridFunction & gf(*_local_test_vars.at(0));
  gf = 0.0;
  mfem::ParLinearForm lf(gf.ParFESpace());
  lf = 0.0;
  mfem::Array<int> ess_bdr_tdofs;
  _bc_map.ApplyEssentialBCs(_h_curl_var_name, ess_bdr_tdofs, gf, _pmesh);
  _bc_map.ApplyIntegratedBCs(_h_curl_var_name, lf, _pmesh);
  lf.Assemble();
  _sources.Apply(&lf);
  mfem::ParBilinearForm blf(gf.ParFESpace());
  blf.AddDomainIntegrator(new mfem::CurlCurlIntegrator(*_stiff_coef));
  blf.Assemble();
  blf.Finalize();
  mfem::HypreParMatrix curl_mu_inv_curl;
  mfem::HypreParVector sol_tdofs(gf.ParFESpace());
  mfem::HypreParVector rhs_tdofs(gf.ParFESpace());
  blf.FormLinearSystem(ess_bdr_tdofs, gf, lf, curl_mu_inv_curl, sol_tdofs, rhs_tdofs);

  // Define and apply a parallel FGMRES solver for AX=B with the AMS
  // preconditioner from hypre.
  hephaestus::DefaultHCurlFGMRESSolver a1_solver(
      _solver_options, curl_mu_inv_curl, gf.ParFESpace());
  a1_solver.Mult(rhs_tdofs, sol_tdofs);
  blf.RecoverFEMSolution(sol_tdofs, lf, gf);
}

} // namespace hephaestus
