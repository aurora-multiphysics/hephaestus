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
    : myid_(0), num_procs_(1), pmesh_(&pmesh), _order(order),
      _fespaces(fespaces), _variables(variables), _bc_map(bc_map),
      _sources(sources), _domain_properties(domain_properties),
      _solver_options(solver_options), a1_solver(NULL) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);

  state_var_names.resize(1);
  state_var_names.at(0) = "h_curl_var";

  aux_var_names.resize(1);
  aux_var_names.at(0) = "curl h_curl_var";
}

void HCurlSolver::Init(mfem::Vector &X) {
  // Define material property coefficients
  SetMaterialCoefficients(_domain_properties);

  _sources.Init(_variables, _fespaces, _bc_map, _domain_properties);

  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    local_test_vars.at(ind)->MakeRef(local_test_vars.at(ind)->ParFESpace(),
                                     const_cast<mfem::Vector &>(X),
                                     true_offsets[ind]);
    *(local_test_vars.at(ind)) = 0.0;
    *(local_trial_vars.at(ind)) = 0.0;
  }

  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("TestVariableName", state_var_names.at(0));
  weak_form_params.SetParam("AlphaCoefName", alpha_coef_name);
  weak_form_params.SetParam("BetaCoefName", beta_coef_name);

  _weak_form = new hephaestus::CurlCurlWeakForm(weak_form_params);
  _weak_form->Init(_variables, _fespaces, _bc_map, _domain_properties);
  _weak_form->buildWeakForm(_bc_map, _sources);
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

  _weak_form->setTimeStep(dt);
  _weak_form->updateWeakForm(_bc_map, _sources);
  _weak_form->FormLinearSystem(A1, X1, B1);
  if (a1_solver == NULL) {
    a1_solver = new hephaestus::DefaultHCurlPCGSolver(_solver_options, A1,
                                                      _weak_form->test_pfes);
  }
  a1_solver->Mult(B1, X1);
  _weak_form->RecoverFEMSolution(X1, *local_trial_vars.at(0));
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

void HCurlSolver::RegisterVariables() {
  RegisterMissingVariables();
  local_test_vars = populateVectorFromNamedFieldsMap<mfem::ParGridFunction>(
      _variables, state_var_names);
  local_trial_vars = registerTimeDerivatives(state_var_names, _variables);

  // Set operator size and block structure
  true_offsets.SetSize(local_test_vars.size() + 1);
  true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    true_offsets[ind + 1] = local_test_vars.at(ind)->ParFESpace()->GetVSize();
  }
  true_offsets.PartialSum();

  this->height = true_offsets[local_test_vars.size()];
  this->width = true_offsets[local_test_vars.size()];

  HYPRE_BigInt state_dofs = true_offsets[local_test_vars.size()];
  if (myid_ == 0) {
    std::cout << "Total number of         DOFs: " << state_dofs << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "Total number of H(Curl) DOFs: " << state_dofs << std::endl;
    std::cout << "------------------------------------" << std::endl;
  }

  // Populate vector of active auxiliary variables
  active_aux_var_names.resize(0);
  for (auto &aux_var_name : aux_var_names) {
    if (_variables.Has(aux_var_name)) {
      active_aux_var_names.push_back(aux_var_name);
    }
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

} // namespace hephaestus
