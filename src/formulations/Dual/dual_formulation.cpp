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

#include "dual_formulation.hpp"

namespace hephaestus {

DualFormulation::DualFormulation() : TimeDomainFormulation() {
  alpha_coef_name = std::string("alpha");
  beta_coef_name = std::string("beta");
  h_curl_var_name = std::string("h_curl_var");
  h_div_var_name = std::string("h_div_var");
}

void DualFormulation::ConstructEquationSystem() {
  hephaestus::InputParameters weak_form_params;
  weak_form_params.SetParam("HCurlVarName", h_curl_var_name);
  weak_form_params.SetParam("HDivVarName", h_div_var_name);
  weak_form_params.SetParam("AlphaCoefName", alpha_coef_name);
  weak_form_params.SetParam("BetaCoefName", beta_coef_name);
  this->GetProblem()->td_equation_system =
      std::make_unique<hephaestus::WeakCurlEquationSystem>(weak_form_params);
}

void DualFormulation::ConstructOperator() {
  this->problem->td_operator = std::make_unique<hephaestus::DualOperator>(
      *(this->problem->pmesh), this->problem->fespaces,
      this->problem->gridfunctions, this->problem->bc_map,
      this->problem->domain_properties, this->problem->sources,
      this->problem->solver_options);
  this->problem->td_operator->SetEquationSystem(
      this->problem->td_equation_system.get());
  this->problem->td_operator->SetVariables();
};

void DualFormulation::RegisterGridFunctions() {
  int &myid = this->GetProblem()->myid_;
  hephaestus::GridFunctions &variables = this->GetProblem()->gridfunctions;

  // Register default ParGridFunctions of state variables if not provided
  if (!variables.Has(h_curl_var_name)) {
    if (myid == 0) {
      MFEM_WARNING(
          h_curl_var_name
          << " not found in variables: building gridfunction from defaults");
    }
    AddFESpace(std::string("_HCurlFESpace"), std::string("ND_3D_P2"));
    AddGridFunction(h_curl_var_name, std::string("_HCurlFESpace"));
  }

  // Register default ParGridFunctions of state variables if not provided
  if (!variables.Has(h_div_var_name)) {
    if (myid == 0) {
      MFEM_WARNING(
          h_div_var_name
          << " not found in variables: building gridfunction from defaults");
    }
    AddFESpace(std::string("_HDivFESpace"), std::string("RT_3D_P3"));
    AddGridFunction(h_div_var_name, std::string("_HDivFESpace"));
  }

  // Register time derivatives
  TimeDomainProblemBuilder::RegisterGridFunctions();
};

void DualFormulation::RegisterCoefficients() {
  hephaestus::DomainProperties &domain_properties =
      this->GetProblem()->domain_properties;
  if (domain_properties.scalar_property_map.count("alpha") == 0) {
    domain_properties.scalar_property_map["alpha"] = new mfem::PWCoefficient(
        domain_properties.getGlobalScalarProperty(std::string("alpha")));
  }
  if (domain_properties.scalar_property_map.count("beta") == 0) {
    domain_properties.scalar_property_map["beta"] = new mfem::PWCoefficient(
        domain_properties.getGlobalScalarProperty(std::string("beta")));
  }
}

WeakCurlEquationSystem::WeakCurlEquationSystem(
    const hephaestus::InputParameters &params)
    : TimeDependentEquationSystem(params),
      h_curl_var_name(params.GetParam<std::string>("HCurlVarName")),
      h_div_var_name(params.GetParam<std::string>("HDivVarName")),
      alpha_coef_name(params.GetParam<std::string>("AlphaCoefName")),
      beta_coef_name(params.GetParam<std::string>("BetaCoefName")),
      dtalpha_coef_name(std::string("dt_") + alpha_coef_name) {}

void WeakCurlEquationSystem::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {
  domain_properties.scalar_property_map[dtalpha_coef_name] =
      new mfem::TransformedCoefficient(
          &dtCoef, domain_properties.scalar_property_map[alpha_coef_name],
          prodFunc);
  TimeDependentEquationSystem::Init(variables, fespaces, bc_map,
                                    domain_properties);
}

void WeakCurlEquationSystem::addKernels() {
  addVariableNameIfMissing(h_curl_var_name);
  addVariableNameIfMissing(h_div_var_name);
  std::string dh_curl_var_dt = GetTimeDerivativeName(h_curl_var_name);

  // (αv_{n}, ∇×u')
  hephaestus::InputParameters weakCurlParams;
  weakCurlParams.SetParam("HCurlVarName", h_curl_var_name);
  weakCurlParams.SetParam("HDivVarName", h_div_var_name);
  weakCurlParams.SetParam("CoefficientName", alpha_coef_name);
  addKernel(h_curl_var_name, new hephaestus::WeakCurlKernel(weakCurlParams));

  // (αdt∇×u_{n+1}, ∇×u')
  hephaestus::InputParameters curlCurlParams;
  curlCurlParams.SetParam("CoefficientName", dtalpha_coef_name);
  addKernel(h_curl_var_name, new hephaestus::CurlCurlKernel(curlCurlParams));

  // (βu_{n+1}, u')
  hephaestus::InputParameters vectorFEMassParams;
  vectorFEMassParams.SetParam("CoefficientName", beta_coef_name);
  addKernel(h_curl_var_name,
            new hephaestus::VectorFEMassKernel(vectorFEMassParams));
}

DualOperator::DualOperator(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : TimeDomainEquationSystemOperator(pmesh, fespaces, variables, bc_map,
                                       domain_properties, sources,
                                       solver_options) {}

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

void DualOperator::SetVariables() {
  TimeDomainEquationSystemOperator::SetVariables();
  // Blocks for solution vector are smaller than the operator size
  // for DualOperator, as curl is stored separately.
  // Block operator only has the HCurl TrueVSize;
  block_trueOffsets.SetSize(local_test_vars.size());
  block_trueOffsets[0] = 0;
  for (unsigned int ind = 0; ind < local_test_vars.size() - 1; ++ind) {
    block_trueOffsets[ind + 1] =
        local_test_vars.at(ind)->ParFESpace()->TrueVSize();
  }
  block_trueOffsets.PartialSum();

  trueX.Update(block_trueOffsets);
  trueRhs.Update(block_trueOffsets);
}

} // namespace hephaestus