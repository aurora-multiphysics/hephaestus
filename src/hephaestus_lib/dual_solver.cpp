// Solves the equations
// ∇⋅s0 = 0
// ∇×(αv) - βu = s0
// dv/dt = -∇×u

// where
// s0 ∈ H(div) source field
// v ∈ H(div)
// u ∈ H(curl)
// p ∈ H1

// Weak form (Space discretisation)
// -(s0, ∇ p') + <n.s0, p'> = 0
// (αv, ∇×u') - (βu, u') - (s0, u') - <(αv) × n, u'> = 0
// (dv/dt, v') + (∇×u, v') = 0

// Time discretisation using implicit scheme:
// Unknowns
// s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
// dv/dt_{n+1} ∈ H(div)
// u_{n+1} ∈ H(curl)
// p_{n+1} ∈ H1

// Fully discretised equations
// -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
// (αv_{n}, ∇×u') - (αdt∇×u_{n+1}, ∇×u') - (βu_{n+1}, u') - (s0_{n+1}, u') -
// <(αv) × n, u'> = 0
// (dv/dt_{n+1}, v') + (∇×u_{n+1}, v') = 0
// using
// v_{n+1} = v_{n} + dt dv/dt_{n+1} = v_{n} - dt ∇×u_{n+1}

// Rewritten as
// a0(p_{n+1}, p') = b0(p')
// a1(u_{n+1}, u') = b1(u')
// dv/dt_{n+1} = -∇×u

// where
// a0(p, p') = (β ∇ p, ∇ p')
// b0(p') = <n.s0, p'>
// a1(u, u') = (βu, u') + (αdt∇×u, ∇×u')
// b1(u') = (s0_{n+1}, u') + (αv_{n}, ∇×u') + <(αdt∇×u_{n+1}) × n, u'>

#include "dual_solver.hpp"

namespace hephaestus {

DualFormulation::DualFormulation() : TransientFormulation() {
  alpha_coef_name = std::string("alpha");
  beta_coef_name = std::string("beta");
  h_curl_var_name = std::string("h_curl_var");
  h_div_var_name = std::string("h_div_var");
}

hephaestus::TimeDependentEquationSystem *
DualFormulation::CreateEquationSystem() {
  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("HCurlVarName", h_curl_var_name);
  weak_form_params.SetParam("HDivVarName", h_div_var_name);
  weak_form_params.SetParam("AlphaCoefName", alpha_coef_name);
  weak_form_params.SetParam("BetaCoefName", beta_coef_name);
  equation_system = new hephaestus::WeakCurlEquationSystem(weak_form_params);
  return equation_system;
}

hephaestus::TimeDomainEquationSystemOperator *
DualFormulation::CreateTimeDomainOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options) {
  td_operator =
      new hephaestus::DualOperator(pmesh, order, fespaces, variables, bc_map,
                                   domain_properties, sources, solver_options);
  td_operator->SetEquationSystem(equation_system);

  return td_operator;
};

void DualFormulation::RegisterMissingVariables(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) {
  int myid;
  MPI_Comm_rank(pmesh.GetComm(), &myid);

  // Register default ParGridFunctions of state variables if not provided
  if (!variables.Has(h_curl_var_name)) {
    if (myid == 0) {
      std::cout
          << h_curl_var_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    fespaces.Register(
        "_HCurlFESpace",
        new mfem::common::ND_ParFESpace(&pmesh, 2, pmesh.Dimension()), true);
    variables.Register(h_curl_var_name,
                       new mfem::ParGridFunction(fespaces.Get("_HCurlFESpace")),
                       true);
  };
  variables.Register(
      hephaestus::TimeDependentEquationSystem::GetTimeDerivativeName(
          h_curl_var_name),
      new mfem::ParGridFunction(fespaces.Get("_HCurlFESpace")), true);

  // Register default ParGridFunctions of state variables if not provided
  if (!variables.Has(h_div_var_name)) {
    if (myid == 0) {
      std::cout
          << h_div_var_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    fespaces.Register(
        "_HDivFESpace",
        new mfem::common::RT_ParFESpace(&pmesh, 2, pmesh.Dimension()), true);
    variables.Register(h_div_var_name,
                       new mfem::ParGridFunction(fespaces.Get("_HDivFESpace")),
                       true);
  };
  variables.Register(
      hephaestus::TimeDependentEquationSystem::GetTimeDerivativeName(
          h_div_var_name),
      new mfem::ParGridFunction(fespaces.Get("_HDivFESpace")), true);
};

void DualFormulation::RegisterCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  if (domain_properties.scalar_property_map.count("alpha") == 0) {
    domain_properties.scalar_property_map["alpha"] = new mfem::PWCoefficient(
        domain_properties.getGlobalScalarProperty(std::string("alpha")));
  }
  if (domain_properties.scalar_property_map.count("beta") == 0) {
    domain_properties.scalar_property_map["beta"] = new mfem::PWCoefficient(
        domain_properties.getGlobalScalarProperty(std::string("beta")));
  }
}

void DualOperator::SetVariables() {
  state_var_names = _equation_system->var_names;
  local_test_vars = populateVectorFromNamedFieldsMap<mfem::ParGridFunction>(
      _variables, _equation_system->var_names);
  local_trial_vars = populateVectorFromNamedFieldsMap<mfem::ParGridFunction>(
      _variables, _equation_system->var_time_derivative_names);

  // Set operator size and block structure
  block_trueOffsets.SetSize(local_test_vars.size());
  block_trueOffsets[0] = 0;
  for (unsigned int ind = 0; ind < local_test_vars.size() - 1; ++ind) {
    block_trueOffsets[ind + 1] =
        local_test_vars.at(ind)->ParFESpace()->TrueVSize();
  }
  block_trueOffsets.PartialSum();

  true_offsets.SetSize(local_test_vars.size() + 1);
  true_offsets[0] = 0;
  for (unsigned int ind = 0; ind < local_test_vars.size(); ++ind) {
    true_offsets[ind + 1] = local_test_vars.at(ind)->ParFESpace()->GetVSize();
  }
  true_offsets.PartialSum();

  this->height = true_offsets[local_test_vars.size()];
  this->width = true_offsets[local_test_vars.size()];
  trueX.Update(block_trueOffsets);
  trueRhs.Update(block_trueOffsets);

  // Populate vector of active auxiliary variables
  active_aux_var_names.resize(0);
  for (auto &aux_var_name : aux_var_names) {
    if (_variables.Has(aux_var_name)) {
      active_aux_var_names.push_back(aux_var_name);
    }
  }
};

void DualOperator::Init(mfem::Vector &X) {
  TimeDomainEquationSystemOperator::Init(X);
  hephaestus::WeakCurlEquationSystem *eqs =
      dynamic_cast<hephaestus::WeakCurlEquationSystem *>(_equation_system);
  h_curl_var_name = eqs->h_curl_var_name;
  h_div_var_name = eqs->h_div_var_name;
  u_ = _variables.Get(h_curl_var_name);
  dv_ = _variables.Get(GetTimeDerivativeName(h_div_var_name));
  HCurlFESpace_ = u_->ParFESpace();
  HDivFESpace_ = dv_->ParFESpace();
  curl = new mfem::ParDiscreteLinearOperator(HCurlFESpace_, HDivFESpace_);
  curl->AddDomainInterpolator(new mfem::CurlInterpolator);
  curl->Assemble();
}

DualOperator::DualOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : TimeDomainEquationSystemOperator(pmesh, order, fespaces, variables,
                                       bc_map, domain_properties, sources,
                                       solver_options) {}

void DualOperator::ImplicitSolve(const double dt, const mfem::Vector &X,
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

  // Subtract off contribution from source
  _sources.SubtractSources(u_);

  // dv/dt_{n+1} = -∇×u
  curl->Mult(*u_, *dv_);
  *dv_ *= -1.0;
}

} // namespace hephaestus
