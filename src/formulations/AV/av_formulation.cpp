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

AVFormulation::AVFormulation(const std::string &alpha_coef_name,
                             const std::string &inv_alpha_coef_name,
                             const std::string &beta_coef_name,
                             const std::string &vector_potential_name,
                             const std::string &scalar_potential_name)
    : TimeDomainEMFormulation(), _alpha_coef_name(alpha_coef_name),
      _inv_alpha_coef_name(inv_alpha_coef_name),
      _beta_coef_name(beta_coef_name),
      _vector_potential_name(vector_potential_name),
      _scalar_potential_name(scalar_potential_name) {}

void AVFormulation::ConstructEquationSystem() {
  hephaestus::InputParameters av_system_params;
  av_system_params.SetParam("VectorPotentialName", _vector_potential_name);
  av_system_params.SetParam("ScalarPotentialName", _scalar_potential_name);
  av_system_params.SetParam("AlphaCoefName", _alpha_coef_name);
  av_system_params.SetParam("BetaCoefName", _beta_coef_name);
  GetProblem()->td_equation_system =
      std::make_unique<hephaestus::AVEquationSystem>(av_system_params);
}

void AVFormulation::RegisterGridFunctions() {
  int &myid = GetProblem()->myid_;
  hephaestus::GridFunctions &gridfunctions = GetProblem()->gridfunctions;

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(_vector_potential_name)) {
    if (myid == 0) {
      MFEM_WARNING(_vector_potential_name
                   << " not found in gridfunctions: building gridfunction from "
                      "defaults");
    }
    AddFESpace(std::string("_HCurlFESpace"), std::string("ND_3D_P2"));
    AddGridFunction(_vector_potential_name, std::string("_HCurlFESpace"));
  }

  // Register default ParGridFunctions of state gridfunctions if not provided
  if (!gridfunctions.Has(_scalar_potential_name)) {
    if (myid == 0) {
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

void AVFormulation::RegisterCoefficients() {
  hephaestus::Coefficients &coefficients = GetProblem()->coefficients;
  if (!coefficients.scalars.Has(_inv_alpha_coef_name)) {
    MFEM_ABORT(_inv_alpha_coef_name + " coefficient not found.");
  }
  if (!coefficients.scalars.Has(_beta_coef_name)) {
    MFEM_ABORT(_beta_coef_name + " coefficient not found.");
  }

  coefficients.scalars.Register(
      _alpha_coef_name,
      new mfem::TransformedCoefficient(
          &oneCoef, coefficients.scalars.Get(_inv_alpha_coef_name), fracFunc),
      true);
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

void AVEquationSystem::Init(hephaestus::GridFunctions &gridfunctions,
                            const hephaestus::FESpaces &fespaces,
                            hephaestus::BCMap &bc_map,
                            hephaestus::Coefficients &coefficients) {
  coefficients.scalars.Register(
      dtalpha_coef_name,
      new mfem::TransformedCoefficient(
          &dtCoef, coefficients.scalars.Get(alpha_coef_name), prodFunc),
      true);

  coefficients.scalars.Register(
      neg_beta_coef_name,
      new mfem::TransformedCoefficient(
          &negCoef, coefficients.scalars.Get(beta_coef_name), prodFunc),
      true);

  TimeDependentEquationSystem::Init(gridfunctions, fespaces, bc_map,
                                    coefficients);
}

void AVEquationSystem::addKernels() {
  addTrialVariableNameIfMissing(a_name);
  std::string da_dt_name = GetTimeDerivativeName(a_name);
  addTrialVariableNameIfMissing(v_name);
  std::string dv_dt_name = GetTimeDerivativeName(v_name);

  // (α∇×A_{n}, ∇×dA'/dt)
  hephaestus::InputParameters weakCurlCurlParams;
  weakCurlCurlParams.SetParam("CoupledVariableName", a_name);
  weakCurlCurlParams.SetParam("CoefficientName", alpha_coef_name);
  addKernel(da_dt_name, new hephaestus::WeakCurlCurlKernel(weakCurlCurlParams));

  // (αdt∇×dA/dt_{n+1}, ∇×dA'/dt)
  hephaestus::InputParameters curlCurlParams;
  curlCurlParams.SetParam("CoefficientName", dtalpha_coef_name);
  addKernel(da_dt_name, new hephaestus::CurlCurlKernel(curlCurlParams));

  // (βdA/dt_{n+1}, dA'/dt)
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

} // namespace hephaestus
