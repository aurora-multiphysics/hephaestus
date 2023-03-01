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

#include "av_solver.hpp"

namespace hephaestus {

AVFormulation::AVFormulation() : TransientFormulation() {
  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");
  CreateEquationSystem();
}

hephaestus::TimeDependentEquationSystem *AVFormulation::CreateEquationSystem() {
  std::vector<std::string> state_var_names;
  state_var_names.resize(2);
  state_var_names.at(0) = "magnetic_vector_potential";
  state_var_names.at(1) = "electric_potential";

  hephaestus::InputParameters av_system_params;
  av_system_params.SetParam("VariableNames", state_var_names);
  av_system_params.SetParam("TimeDerivativeNames",
                            GetTimeDerivativeNames(state_var_names));
  av_system_params.SetParam("AlphaCoefName", alpha_coef_name);
  av_system_params.SetParam("BetaCoefName", beta_coef_name);
  equation_system = new hephaestus::AVEquationSystem(av_system_params);
  return equation_system;
}

hephaestus::TimeDomainEquationSystemOperator *
AVFormulation::CreateTimeDomainOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options) {
  td_operator =
      new hephaestus::AVOperator(pmesh, order, fespaces, variables, bc_map,
                                 domain_properties, sources, solver_options);
  td_operator->SetEquationSystem(equation_system);

  return td_operator;
};

void AVFormulation::RegisterMissingVariables(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) {
  int myid;
  MPI_Comm_rank(pmesh.GetComm(), &myid);

  // Register default ParGridFunctions of state variables if not provided
  std::string u_name = equation_system->var_names.at(0);
  if (!variables.Has(u_name)) {
    if (myid == 0) {
      std::cout
          << u_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    fespaces.Register(
        "_HCurlFESpace",
        new mfem::common::ND_ParFESpace(&pmesh, 2, pmesh.Dimension()), true);
    variables.Register(
        u_name, new mfem::ParGridFunction(fespaces.Get("_HCurlFESpace")), true);
  };
  std::string p_name = equation_system->var_names.at(1);
  if (!variables.Has(p_name)) {
    if (myid == 0) {
      std::cout
          << p_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    fespaces.Register(
        "_H1FESpace",
        new mfem::common::H1_ParFESpace(&pmesh, 2, pmesh.Dimension()), true);
    variables.Register(
        p_name, new mfem::ParGridFunction(fespaces.Get("_H1FESpace")), true);
  };

  RegisterTimeDerivatives(equation_system->var_names, variables);
};

void AVFormulation::RegisterCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  if (domain_properties.scalar_property_map.count("magnetic_permeability") ==
      0) {
    domain_properties.scalar_property_map["magnetic_permeability"] =
        new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
            std::string("magnetic_permeability")));
  }
  if (domain_properties.scalar_property_map.count("electrical_conductivity") ==
      0) {
    domain_properties.scalar_property_map["electrical_conductivity"] =
        new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
            std::string("electrical_conductivity")));
  }

  domain_properties.scalar_property_map[alpha_coef_name] =
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map["magnetic_permeability"],
          fracFunc);
}

AVOperator::AVOperator(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : TimeDomainEquationSystemOperator(pmesh, order, fespaces, variables,
                                       bc_map, domain_properties, sources,
                                       solver_options) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);

  aux_var_names.resize(2);
  aux_var_names.at(0) = "electric_field";
  aux_var_names.at(1) = "magnetic_flux_density";
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
void AVOperator::ImplicitSolve(const double dt, const mfem::Vector &X,
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
  if (solver != NULL) {
    delete solver;
  }
  solver = new hephaestus::DefaultGMRESSolver(
      _solver_options, *blockA.As<mfem::HypreParMatrix>());
  // solver = new hephaestus::DefaultGMRESSolver(_solver_options, *blockA,
  //                                             pmesh_->GetComm());

  solver->Mult(trueRhs, trueX);
  _equation_system->RecoverFEMSolution(trueX, _variables);

} // namespace hephaestus
} // namespace hephaestus
