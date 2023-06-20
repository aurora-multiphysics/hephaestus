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

void AVFormulation::ConstructEquationSystem() {
  hephaestus::InputParameters av_system_params;
  av_system_params.SetParam("VectorPotentialName", vector_potential_name);
  av_system_params.SetParam("ScalarPotentialName", scalar_potential_name);
  av_system_params.SetParam("AlphaCoefName", alpha_coef_name);
  av_system_params.SetParam("BetaCoefName", beta_coef_name);
  this->GetProblem()->td_equation_system =
      std::make_unique<hephaestus::AVEquationSystem>(av_system_params);
}

void AVFormulation::ConstructOperator() {
  this->problem->td_operator = std::make_unique<hephaestus::AVOperator>(
      *(this->problem->pmesh), this->problem->fespaces,
      this->problem->gridfunctions, this->problem->bc_map,
      this->problem->domain_properties, this->problem->sources,
      this->problem->solver_options);
  this->problem->td_operator->SetEquationSystem(
      this->problem->td_equation_system.get());
  this->problem->td_operator->SetVariables();
};

void AVFormulation::RegisterGridFunctions() {
  int &myid = this->GetProblem()->myid_;
  hephaestus::GridFunctions &variables = this->GetProblem()->gridfunctions;

  // Register default ParGridFunctions of state variables if not provided
  if (!variables.Has(vector_potential_name)) {
    if (myid == 0) {
      MFEM_WARNING(
          vector_potential_name
          << " not found in variables: building gridfunction from defaults");
    }
    AddFESpace(std::string("_HCurlFESpace"), std::string("ND_3D_P2"));
    AddGridFunction(vector_potential_name, std::string("_HCurlFESpace"));
  }

  // Register default ParGridFunctions of state variables if not provided
  if (!variables.Has(scalar_potential_name)) {
    if (myid == 0) {
      MFEM_WARNING(
          scalar_potential_name
          << " not found in variables: building gridfunction from defaults");
    }
    AddFESpace(std::string("_H1FESpace"), std::string("H1_3D_P2"));
    AddGridFunction(scalar_potential_name, std::string("_H1FESpace"));
  }

  // Register time derivatives
  TimeDomainProblemBuilder::RegisterGridFunctions();
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

AVEquationSystem::AVEquationSystem(const hephaestus::InputParameters &params)
    : TimeDependentEquationSystem(params),
      a_name(params.GetParam<std::string>("VectorPotentialName")),
      v_name(params.GetParam<std::string>("ScalarPotentialName")),
      alpha_coef_name(params.GetParam<std::string>("AlphaCoefName")),
      beta_coef_name(params.GetParam<std::string>("BetaCoefName")),
      dtalpha_coef_name(std::string("dt_") + alpha_coef_name),
      neg_beta_coef_name(std::string("negative_") + beta_coef_name),
      negCoef(-1.0) {}

void AVEquationSystem::Init(
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    const mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    hephaestus::BCMap &bc_map,
    hephaestus::DomainProperties &domain_properties) {
  domain_properties.scalar_property_map[dtalpha_coef_name] =
      new mfem::TransformedCoefficient(
          &dtCoef, domain_properties.scalar_property_map[alpha_coef_name],
          prodFunc);
  domain_properties.scalar_property_map[neg_beta_coef_name] =
      new mfem::TransformedCoefficient(
          &negCoef, domain_properties.scalar_property_map[beta_coef_name],
          prodFunc);
  TimeDependentEquationSystem::Init(variables, fespaces, bc_map,
                                    domain_properties);
}

void AVEquationSystem::addKernels() {
  addVariableNameIfMissing(a_name);
  std::string da_dt_name = GetTimeDerivativeName(a_name);
  addVariableNameIfMissing(v_name);
  std::string dv_dt_name = GetTimeDerivativeName(v_name);

  // (α∇×A_{n}, ∇×A')
  hephaestus::InputParameters weakCurlCurlParams;
  weakCurlCurlParams.SetParam("CoupledVariableName", a_name);
  weakCurlCurlParams.SetParam("CoefficientName", alpha_coef_name);
  addKernel(da_dt_name, new hephaestus::WeakCurlCurlKernel(weakCurlCurlParams));

  // (αdt∇×dA/dt_{n+1}, ∇×A')
  hephaestus::InputParameters curlCurlParams;
  curlCurlParams.SetParam("CoefficientName", dtalpha_coef_name);
  addKernel(da_dt_name, new hephaestus::CurlCurlKernel(curlCurlParams));

  // (βdA/dt_{n+1}, A')
  hephaestus::InputParameters vectorFEMassParams;
  vectorFEMassParams.SetParam("CoefficientName", beta_coef_name);
  addKernel(da_dt_name, new hephaestus::VectorFEMassKernel(vectorFEMassParams));

  // (σ ∇ V, dA'/dt)
  hephaestus::InputParameters mixedVectorGradientParams;
  mixedVectorGradientParams.SetParam("CoefficientName", beta_coef_name);
  addKernel(
      v_name, da_dt_name,
      new hephaestus::MixedVectorGradientKernel(mixedVectorGradientParams));

  // (σ ∇ V, ∇ V')
  hephaestus::InputParameters diffusionParams;
  diffusionParams.SetParam("CoefficientName", beta_coef_name);
  addKernel(v_name, new hephaestus::DiffusionKernel(diffusionParams));

  // (σdA/dt, ∇ V')
  hephaestus::InputParameters vectorFEWeakDivergenceParams;
  vectorFEWeakDivergenceParams.SetParam("CoefficientName", beta_coef_name);
  addKernel(da_dt_name, v_name,
            new hephaestus::VectorFEWeakDivergenceKernel(
                vectorFEWeakDivergenceParams));
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
}
} // namespace hephaestus