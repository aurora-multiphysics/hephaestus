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

AVSolver::AVSolver(mfem::ParMesh &pmesh, int order,
                   mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
                   mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                   hephaestus::BCMap &bc_map,
                   hephaestus::DomainProperties &domain_properties,
                   hephaestus::Sources &sources,
                   hephaestus::InputParameters &solver_options)
    : TransientFormulation(pmesh, order, fespaces, variables, bc_map,
                           domain_properties, sources, solver_options) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);

  state_var_names.resize(2);
  state_var_names.at(0) = "magnetic_vector_potential";
  state_var_names.at(1) = "electric_potential";

  aux_var_names.resize(2);
  aux_var_names.at(0) = "electric_field";
  aux_var_names.at(1) = "magnetic_flux_density";
}

void AVSolver::RegisterMissingVariables() {
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

  std::string p_name = state_var_names.at(1);
  if (!_variables.Has(p_name)) {
    if (myid_ == 0) {
      std::cout
          << p_name
          << " not found in variables: building gridfunction from defaults"
          << std::endl;
    }
    _fespaces.Register(
        "_H1FESpace",
        new mfem::common::H1_ParFESpace(pmesh_, _order, pmesh_->Dimension()),
        true);
    _variables.Register(
        p_name, new mfem::ParGridFunction(_fespaces.Get("_H1FESpace")), true);
  }
}

void AVSolver::SetMaterialCoefficients(
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

  alpha_coef_name = std::string("magnetic_reluctivity");
  beta_coef_name = std::string("electrical_conductivity");

  domain_properties.scalar_property_map[alpha_coef_name] =
      new mfem::TransformedCoefficient(
          &oneCoef,
          domain_properties.scalar_property_map["magnetic_permeability"],
          fracFunc);
}

void AVSolver::SetEquationSystem() {
  hephaestus::InputParameters av_system_params;
  av_system_params.SetParam("VariableNames", state_var_names);
  av_system_params.SetParam("TimeDerivativeNames",
                            GetTimeDerivativeNames(state_var_names));
  av_system_params.SetParam("AlphaCoefName", alpha_coef_name);
  av_system_params.SetParam("BetaCoefName", beta_coef_name);

  _equation_system = new hephaestus::AVEquationSystem(av_system_params);
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
void AVSolver::ImplicitSolve(const double dt, const mfem::Vector &X,
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
  // _equation_system->RecoverFEMSolution(trueX, *local_trial_vars.at(0));

  trueX.GetBlock(0).SyncAliasMemory(trueX);
  trueX.GetBlock(1).SyncAliasMemory(trueX);

  local_trial_vars.at(0)->Distribute(&(trueX.GetBlock(0)));
  local_test_vars.at(1)->Distribute(&(trueX.GetBlock(1)));

  // trueX.GetBlock(0).SyncAliasMemory(trueX);
  // trueX.GetBlock(1).SyncAliasMemory(trueX);

  // du_.Distribute(&(trueX.GetBlock(0)));
  // p_.Distribute(&(trueX.GetBlock(1)));
  // dX_dt = 0.0;
  // dtCoef.constant = dt;

  // u_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  // p_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[1]);

  // du_.MakeRef(HCurlFESpace_, dX_dt, true_offsets[0]);
  // dp_.MakeRef(H1FESpace_, dX_dt, true_offsets[1]);

  // _domain_properties.SetTime(this->GetTime());

  //////////////////////////////////////////////////////////////////////////////
  // mfem::BlockVector trueX(block_trueOffsets), trueRhs(block_trueOffsets);
  // trueX = 0.0;
  // *b0 = 0.0;
  // *b01 = 0.0;
  // *b10 = 0.0;

  // // -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
  // // a0(p_{n+1}, p') = b0(p')
  // // a0(p, p') = (β ∇ p, ∇ p')
  // // b0(p') = <n.s0, p'>

  // (α∇×u_{n}, ∇×u') + (αdt∇×du/dt_{n+1}, ∇×u') + (βdu/dt_{n+1}, u')
  // - (s0_{n+1}, u') - <(α∇×u_{n+1}) × n, u'> = 0
  // a1(du/dt_{n+1}, u') = b1(u')
  // a1(u, u') = (βu, u') + (αdt∇×u, ∇×u')
  // b1(u') = (s0_{n+1}, u') - (α∇×u_{n}, ∇×u') + <(α∇×u_{n+1}) × n, u'>

  // (α∇×u_{n}, ∇×u')
  // curlCurl->MultTranspose(u_, *b0);
  // *b0 *= -1.0;

  // mfem::ParGridFunction J_gf(HCurlFESpace_);
  // mfem::Array<int> ess_tdof_list;
  // J_gf = 0.0;
  // _bc_map.applyEssentialBCs(u_name, ess_tdof_list, J_gf, pmesh_);
  // _bc_map.applyIntegratedBCs(u_name, *b0, pmesh_);
  // b0->Assemble();

  // _sources.Apply(b0);

  // mfem::ParGridFunction Phi_gf(H1FESpace_);
  // mfem::Array<int> poisson_ess_tdof_list;
  // Phi_gf = 0.0;
  // *b1 = 0.0;
  // _bc_map.applyEssentialBCs(p_name, poisson_ess_tdof_list, Phi_gf, pmesh_);
  // _bc_map.applyIntegratedBCs(p_name, *b1, pmesh_);
  // b1->Assemble();

  // if (a0 == NULL || a10 == NULL || fabs(dt - dt_A0) > 1.0e-12 * dt) {
  //   this->buildA0(betaCoef, dtAlphaCoef);
  // }

  // if (a1 == NULL || a01 == NULL || fabs(dt - dt_A1) > 1.0e-12 * dt) {
  //   this->buildA1(betaCoef);
  // }

  // a0->FormLinearSystem(ess_tdof_list, J_gf, *b0, *A0, *X0, *B0);
  // trueX.GetBlock(0) = *X0;
  // trueRhs.GetBlock(0) = *B0;
  // a1->FormLinearSystem(poisson_ess_tdof_list, Phi_gf, *b1, *A1, *X1, *B1);
  // trueX.GetBlock(1) = *X1;
  // trueRhs.GetBlock(1) = *B1;
  // // b01 = new mfem::ParLinearForm(HCurlFESpace_);
  // // b10 = new mfem::ParLinearForm(H1FESpace_);

  // a01->FormRectangularLinearSystem(ess_tdof_list, poisson_ess_tdof_list,
  // J_gf,
  //                                  *b10, *A01, *X0, *B1);
  // a10->FormRectangularLinearSystem(poisson_ess_tdof_list, ess_tdof_list,
  // Phi_gf,
  //                                  *b01, *A10, *X1, *B0);
  // trueRhs.GetBlock(1) += *B1;
  // trueRhs.GetBlock(0) += *B0;

  // trueX.GetBlock(0).SyncAliasMemory(trueX);
  // trueX.GetBlock(1).SyncAliasMemory(trueX);

  // trueRhs.GetBlock(0).SyncAliasMemory(trueRhs);
  // trueRhs.GetBlock(1).SyncAliasMemory(trueRhs);

  // mfem::Array2D<mfem::HypreParMatrix *> hBlocks(2, 2);
  // hBlocks = NULL;
  // hBlocks(0, 0) = A0;
  // hBlocks(0, 1) = A10;
  // hBlocks(1, 0) = A01;
  // hBlocks(1, 1) = A1;
  // blockA = mfem::HypreParMatrixFromBlocks(hBlocks);

  // if (solver == NULL) {
  //   solver = new hephaestus::DefaultGMRESSolver(_solver_options, *blockA,
  //                                               pmesh_->GetComm());
  // }
  // solver->Mult(trueRhs, trueX);

  // #ifdef MFEM_USE_MUMPS
  //   mfem::MUMPSSolver mumps;
  //   mumps.SetPrintLevel(0);
  //   mumps.SetMatrixSymType(mfem::MUMPSSolver::MatType::UNSYMMETRIC);
  //   mumps.SetOperator(*blockA);
  //   mumps.Mult(trueRhs, trueX);
  // #else
  //   mfem::GMRESSolver gmres(HCurlFESpace_->GetComm());
  //   gmres.SetOperator(*blockA);
  //   gmres.SetAbsTol(1e-16);
  //   gmres.SetMaxIter(10000);
  //   gmres.Mult(trueRhs, trueX);
  // #endif

  // trueX.GetBlock(0).SyncAliasMemory(trueX);
  // trueX.GetBlock(1).SyncAliasMemory(trueX);

  // local_trial_vars.at(0)->Distribute(&(trueX.GetBlock(0)));
  // local_test_vars.at(1)->Distribute(&(trueX.GetBlock(1)));

  // du_.Distribute(&(trueX.GetBlock(0)));
  // p_.Distribute(&(trueX.GetBlock(1)));

  // grad->Mult(p_, e_);
  // e_ += du_;
  // e_ *= -1.0;

  // curl->Mult(u_, b_);
  // curl->AddMult(du_, b_, dt);
}

// void AVSolver::buildA0(mfem::Coefficient *betaCoef,
//                        mfem::Coefficient *dtAlphaCoef) {
//   if (a0 != NULL) {
//     delete a0;
//   }
//   if (a10 != NULL) {
//     delete a10;
//   }

//   // a0(dA/dt, dA'/dt) = (βdA/dt_{n+1}, dA'/dt) + (αdt∇×dA/dt_{n+1},
//   ∇×dA'/dt) a0 = new mfem::ParBilinearForm(HCurlFESpace_);
//   a0->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*betaCoef));
//   a0->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*dtAlphaCoef));
//   a0->Assemble();

//   // a10(V, dA'/dt) = (σ ∇ V, dA'/dt)
//   a10 = new mfem::ParMixedBilinearForm(H1FESpace_, HCurlFESpace_);
//   a10->AddDomainIntegrator(new
//   mfem::MixedVectorGradientIntegrator(*betaCoef)); a10->Assemble();

//   dt_A0 = dtCoef.constant;
// }

// void AVSolver::buildA1(mfem::Coefficient *betaCoef) {
//   if (a1 != NULL) {
//     delete a1;
//   }
//   if (a01 != NULL) {
//     delete a01;
//   }
//   if (negBetaCoef != NULL) {
//     delete negBetaCoef;
//   }

//   negCoef = mfem::ConstantCoefficient(-1.0);
//   negBetaCoef = new mfem::TransformedCoefficient(&negCoef, betaCoef,
//   prodFunc);

//   // a1(V, V') = (σ ∇ V, ∇ V')
//   a1 = new mfem::ParBilinearForm(H1FESpace_);
//   a1->AddDomainIntegrator(new mfem::DiffusionIntegrator(*negBetaCoef));
//   a1->Assemble();

//   // (σdA/dt, ∇ V')
//   a01 = new mfem::ParMixedBilinearForm(HCurlFESpace_, H1FESpace_);
//   a01->AddDomainIntegrator(
//       new mfem::VectorFEWeakDivergenceIntegrator(*negBetaCoef));
//   a01->Assemble();

//   dt_A1 = dtCoef.constant;
// }

// void AVSolver::buildGrad() {
//   // Discrete Grad operator
//   if (grad != NULL) {
//     delete grad;
//   }
//   grad = new mfem::ParDiscreteLinearOperator(H1FESpace_, HCurlFESpace_);
//   grad->AddDomainInterpolator(new mfem::GradientInterpolator());
//   grad->Assemble();
// }

// void AVSolver::buildCurl(mfem::Coefficient *MuInv) {
//   if (curlCurl != NULL) {
//     delete curlCurl;
//   }
//   curlCurl = new mfem::ParBilinearForm(HCurlFESpace_);
//   curlCurl->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*MuInv));
//   curlCurl->Assemble();

//   // Discrete Curl operator
//   if (curl != NULL) {
//     delete curl;
//   }
//   curl = new mfem::ParDiscreteLinearOperator(HCurlFESpace_, HDivFESpace_);
//   curl->AddDomainInterpolator(new mfem::CurlInterpolator());
//   curl->Assemble();
// }

} // namespace hephaestus
