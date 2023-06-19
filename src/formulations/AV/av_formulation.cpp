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

namespace hephaestus {

AVFormulation::AVFormulation() : TimeDomainFormulation() {
  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");
  vector_potential_name = std::string("magnetic_vector_potential");
  scalar_potential_name = std::string("electric_potential");
}

std::unique_ptr<hephaestus::TimeDependentEquationSystem>
AVFormulation::CreateTimeDependentEquationSystem() const {
  hephaestus::InputParameters av_system_params;
  av_system_params.SetParam("VectorPotentialName", vector_potential_name);
  av_system_params.SetParam("ScalarPotentialName", scalar_potential_name);
  av_system_params.SetParam("AlphaCoefName", alpha_coef_name);
  av_system_params.SetParam("BetaCoefName", beta_coef_name);
  return std::make_unique<hephaestus::AVEquationSystem>(av_system_params);
}

std::unique_ptr<hephaestus::TimeDomainEquationSystemOperator>
AVFormulation::CreateTimeDomainEquationSystemOperator(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources,
    hephaestus::InputParameters &solver_options) const {
  return std::make_unique<hephaestus::AVOperator>(pmesh, fespaces, variables,
                                                  bc_map, domain_properties,
                                                  sources, solver_options);
};

void AVFormulation::RegisterMissingVariables(
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables) {
  int myid;
  MPI_Comm_rank(pmesh.GetComm(), &myid);

  // Register default ParGridFunctions of state variables if not provided
  if (!variables.Has(vector_potential_name)) {
    if (myid == 0) {
      std::cout
          << vector_potential_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    fespaces.Register(
        "_HCurlFESpace",
        new mfem::common::ND_ParFESpace(&pmesh, 2, pmesh.Dimension()), true);
    variables.Register(vector_potential_name,
                       new mfem::ParGridFunction(fespaces.Get("_HCurlFESpace")),
                       true);
  };
  variables.Register(
      hephaestus::TimeDependentEquationSystem::GetTimeDerivativeName(
          vector_potential_name),
      new mfem::ParGridFunction(
          variables.Get(vector_potential_name)->ParFESpace()),
      true);
  if (!variables.Has(scalar_potential_name)) {
    if (myid == 0) {
      std::cout
          << scalar_potential_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    fespaces.Register(
        "_H1FESpace",
        new mfem::common::H1_ParFESpace(&pmesh, 2, pmesh.Dimension()), true);
    variables.Register(scalar_potential_name,
                       new mfem::ParGridFunction(fespaces.Get("_H1FESpace")),
                       true);
  };
  variables.Register(
      hephaestus::TimeDependentEquationSystem::GetTimeDerivativeName(
          scalar_potential_name),
      new mfem::ParGridFunction(
          variables.Get(scalar_potential_name)->ParFESpace()),
      true);
};

void AVFormulation::RegisterCoefficients() {
  hephaestus::DomainProperties &domain_properties =
      this->GetProblem()->domain_properties;
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
    mfem::ParMesh &pmesh,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : TimeDomainEquationSystemOperator(pmesh, fespaces, variables, bc_map,
                                       domain_properties, sources,
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
