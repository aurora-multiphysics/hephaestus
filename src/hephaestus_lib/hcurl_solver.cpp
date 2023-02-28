// Solves the equations
// ∇⋅s0 = 0
// ∇×(α∇×u) + βdu/dt = s0

// where
// s0 ∈ H(div) source field
// u ∈ H(curl)
// p ∈ H1

// Dirichlet boundaries constrain du/dt
// Integrated boundaries constrain (α∇×u) × n

// Weak form (Space discretisation)
// -(s0, ∇ p') + <n.s0, p'> = 0
// (α∇×u, ∇×u') + (βdu/dt, u') - (s0, u') - <(α∇×u) × n, u'> = 0

// Time discretisation using implicit scheme:
// Unknowns
// s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
// du/dt_{n+1} ∈ H(curl)
// p_{n+1} ∈ H1

// Fully discretised equations
// -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
// (α∇×u_{n}, ∇×u') + (αdt∇×du/dt_{n+1}, ∇×u') + (βdu/dt_{n+1}, u')
// - (s0_{n+1}, u') - <(α∇×u_{n+1}) × n, u'> = 0
// using
// u_{n+1} = u_{n} + dt du/dt_{n+1}

// Rewritten as
// a0(p_{n+1}, p') = b0(p')
// a1(du/dt_{n+1}, u') = b1(u')

// where
// a0(p, p') = (β ∇ p, ∇ p')
// b0(p') = <n.s0, p'>
// a1(u, u') = (βu, u') + (αdt∇×u, ∇×u')
// b1(u') = (s0_{n+1}, u') - (α∇×u_{n}, ∇×u') + <(α∇×u_{n+1}) × n, u'>

// Formulation should add equation system(s) and required variables
// - Add state equation system
// Formulation - factory class for building:
// -- state variables if not present
// -- auxkernels if auxvariables exist
// -- references to required coefficients
// -- TimeDependentOperator for solve

#include "hcurl_solver.hpp"

namespace hephaestus {

HCurlSolver::HCurlSolver(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : TransientFormulation(pmesh, order, fespaces, variables, bc_map,
                           domain_properties, sources, solver_options) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);

  state_var_names.resize(1);
  state_var_names.at(0) = "h_curl_var";

  aux_var_names.resize(1);
  aux_var_names.at(0) = "curl h_curl_var";
}

/*
This is the main computational code that computes dX/dt implicitly
where X is the state vector containing u.

Unknowns
s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
du/dt_{n+1} ∈ H(curl)
p_{n+1} ∈ H1

Fully discretised equations
(α∇×u_{n}, ∇×u') + (αdt∇×du/dt_{n+1}, ∇×u') + (βdu/dt_{n+1}, u')
- (s0_{n+1}, u') - <(α∇×u_{n+1}) × n, u'> = 0
using
u_{n+1} = u_{n} + dt du/dt_{n+1}
*/
// Should become identical to AVSolver
void HCurlSolver::ImplicitSolve(const double dt, const mfem::Vector &X,
                                mfem::Vector &dX_dt) {
  dX_dt = 0.0;
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    local_test_vars.at(ind)->MakeRef(local_test_vars.at(ind)->ParFESpace(),
                                     const_cast<mfem::Vector &>(X),
                                     true_offsets[ind]);
    local_trial_vars.at(ind)->MakeRef(local_trial_vars.at(ind)->ParFESpace(),
                                      dX_dt, true_offsets[ind]);
  }
  _domain_properties.SetTime(this->GetTime());
  _equation_system->setTimeStep(dt);
  _equation_system->updateEquationSystem(_bc_map, _sources);
  _equation_system->FormLinearSystem(blockA, trueX, trueRhs);
  if (a1_solver != NULL) {
    delete a1_solver;
  }
  a1_solver = new hephaestus::DefaultHCurlPCGSolver(
      _solver_options, *blockA.As<mfem::HypreParMatrix>(),
      _equation_system->test_pfespaces.at(0));
  a1_solver->Mult(trueRhs, trueX);
  _equation_system->RecoverFEMSolution(trueX, _variables);
}

void HCurlSolver::RegisterMissingVariables() {
  // Register default ParGridFunctions of state variables if not provided
  std::string u_name = state_var_names.at(0);
  if (!_variables.Has(u_name)) {
    if (myid_ == 0) {
      std::cout
          << u_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    _fespaces.Register(
        "_HCurlFESpace",
        new mfem::common::ND_ParFESpace(pmesh_, _order, pmesh_->Dimension()),
        true);
    _variables.Register(
        u_name, new mfem::ParGridFunction(_fespaces.Get("_HCurlFESpace")),
        true);
  }
}

void HCurlSolver::SetMaterialCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  if (domain_properties.scalar_property_map.count("alpha") == 0) {
    domain_properties.scalar_property_map["alpha"] = new mfem::PWCoefficient(
        domain_properties.getGlobalScalarProperty(std::string("alpha")));
  }
  if (domain_properties.scalar_property_map.count("beta") == 0) {
    domain_properties.scalar_property_map["beta"] = new mfem::PWCoefficient(
        domain_properties.getGlobalScalarProperty(std::string("beta")));
  }
  alpha_coef_name = std::string("alpha");
  beta_coef_name = std::string("beta");
}

void HCurlSolver::SetEquationSystem() {
  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("VariableNames", state_var_names);
  weak_form_params.SetParam("TimeDerivativeNames",
                            GetTimeDerivativeNames(state_var_names));
  weak_form_params.SetParam("AlphaCoefName", alpha_coef_name);
  weak_form_params.SetParam("BetaCoefName", beta_coef_name);

  _equation_system = new hephaestus::CurlCurlEquationSystem(weak_form_params);
}

} // namespace hephaestus
