// Strong form
// ∇⋅v0 = 0
// ∇×(α∇×u) - βdu/dt = v0

// Weak form
// -(v0, ∇ p') + <n.v0, p'> = 0
// (α∇×u, ∇×u') - (βdu/dt, u') - (v0, u') - <(α∇×u) × n, u'> = 0

// v, v0 ∈ H(div) source field
// u ∈ H(curl)
// p ∈ H1

///////
// (σ(dA/dt + ∇ V), A') +(ν∇×A, ∇×A') - <(ν∇×A) × n, A'> - (J0, A') = 0
// (σ(dA/dt + ∇ V), ∇ V') + <n.J, V'> = 0
///////

// (σ(dA/dt + ∇ V), A') +(ν∇×A, ∇×A') - <(ν∇×A) × n, A'> - (J0, A') = 0

// // (σ dA/dt, A')
// mfem::ParBilinearForm(HCurlFESpace_);
// mfem::VectorFEMassIntegrator(sigmaCoef);

// // (σ ∇ V, A')
// mfem::ParMixedBilinearForm(H1FESpace_, HCurlFESpace_);
// mfem::MixedVectorGradientIntegrator(sigmaCoef);

// // (ν∇×A, ∇×A')
// mfem::ParBilinearForm(HCurlFESpace_);
// mfem::CurlCurlIntegrator(dtMuInvCoef);
// mfem::ParLinearForm(HCurlFESpace_);
// mfem::CurlCurlIntegrator(muInvCoef);

// (σ(dA/dt + ∇ V), ∇ V') + <n.J, V'> = 0

// (σdA/dt , ∇ V')
// mfem::ParMixedBilinearForm(HCurlFESpace_, H1FESpace_)
// mfem::MixedVectorWeakDivergenceIntegrator(sigmaCoef)

// (σ ∇ V, ∇ V')
// mfem::ParBilinearForm(H1FESpace_)
// mfem::DiffusionIntegrator(sigmaCoef)

#include "av_solver.hpp"

namespace hephaestus {

AVSolver::AVSolver(mfem::ParMesh &pmesh, int order,
                   mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                   hephaestus::BCMap &bc_map,
                   hephaestus::DomainProperties &domain_properties)
    : myid_(0), num_procs_(1), pmesh_(&pmesh), _variables(variables),
      _bc_map(bc_map), _domain_properties(domain_properties),
      H1FESpace_(
          new mfem::common::H1_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HCurlFESpace_(
          new mfem::common::ND_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HDivFESpace_(
          new mfem::common::RT_ParFESpace(&pmesh, order, pmesh.Dimension())),
      blockAV(NULL), blockAVPr(NULL), a1(NULL), amg_a0(NULL), pcg_a0(NULL),
      ams_a1(NULL), pcg_a1(NULL), m1(NULL), grad(NULL), curl(NULL),
      weakCurl(NULL), curlCurl(NULL), sourceVecCoef(NULL), src_gf(NULL),
      div_free_src_gf(NULL), hCurlMass(NULL), divFreeProj(NULL),
      v_(mfem::ParGridFunction(H1FESpace_)),
      a_(mfem::ParGridFunction(HCurlFESpace_)),
      e_(mfem::ParGridFunction(HCurlFESpace_)),
      b_(mfem::ParGridFunction(HDivFESpace_)),
      dv_(mfem::ParGridFunction(H1FESpace_)),
      da_(mfem::ParGridFunction(HCurlFESpace_)) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);

  true_offsets.SetSize(3);
  true_offsets[0] = 0;
  true_offsets[1] = H1FESpace_->GetVSize();
  true_offsets[2] = HCurlFESpace_->GetVSize();
  true_offsets.PartialSum();

  block_trueOffsets.SetSize(3);
  block_trueOffsets[0] = 0;
  block_trueOffsets[1] = H1FESpace_->TrueVSize();
  block_trueOffsets[2] = HCurlFESpace_->TrueVSize();
  block_trueOffsets.PartialSum();

  this->height = true_offsets[2];
  this->width = true_offsets[2];
}

void AVSolver::Init(mfem::Vector &X) {
  SetVariableNames();
  _variables.Register(u_name, &a_, false);
  _variables.Register(p_name, &v_, false);
  _variables.Register(v_name, &b_, false);
  _variables.Register(e_name, &e_, false);

  // Define material property coefficients
  dtCoef = mfem::ConstantCoefficient(1.0);
  oneCoef = mfem::ConstantCoefficient(1.0);
  negCoef = mfem::ConstantCoefficient(-1.0);

  SetMaterialCoefficients(_domain_properties);
  dtAlphaCoef = new mfem::TransformedCoefficient(&dtCoef, alphaCoef, prodFunc);
  negBetaCoef = new mfem::TransformedCoefficient(&oneCoef, betaCoef, prodFunc);

  x = new mfem::BlockVector(true_offsets);
  rhs = new mfem::BlockVector(true_offsets);
  trueX = new mfem::BlockVector(block_trueOffsets);
  trueRhs = new mfem::BlockVector(block_trueOffsets);

  // Variables
  // v_, "electric_potential"
  // a_, "magnetic_vector_potential"

  // Coefficients
  // σ, "electrical_conductivity"
  // mu, "magnetic_permeability"

  ///////
  // (σ(dA/dt + ∇ V), ∇ V') + <n.J, V'> = 0
  // (σ(dA/dt + ∇ V), A') + (ν∇×A, ∇×A') - <(ν∇×A) × n, A'> - (J0, A') = 0
  ///////
  // Bilinear for divergence free source field solve
  // -(J0, ∇ V') + <n.J, V'> = 0  where J0 = -σ∇V
  a0 = new mfem::ParBilinearForm(H1FESpace_);
  a0->AddDomainIntegrator(new mfem::DiffusionIntegrator(*betaCoef));
  a0->Assemble();

  // (σdA/dt , ∇ V')
  // mfem::ParMixedBilinearForm(HCurlFESpace_, H1FESpace_)
  // mfem::MixedVectorWeakDivergenceIntegrator(sigmaCoef)
  // TODO: SHOULD BE -BETACOEF
  a01 = new mfem::ParMixedBilinearForm(HCurlFESpace_, H1FESpace_);
  a01->AddDomainIntegrator(
      new mfem::MixedVectorWeakDivergenceIntegrator(*negBetaCoef));
  a01->Assemble();
  // // (σ ∇ V, A')
  // mfem::ParMixedBilinearForm(H1FESpace_, HCurlFESpace_);
  // mfem::MixedVectorGradientIntegrator(sigmaCoef);
  a10 = new mfem::ParMixedBilinearForm(H1FESpace_, HCurlFESpace_);
  a10->AddDomainIntegrator(new mfem::MixedVectorGradientIntegrator(*betaCoef));
  a10->Assemble();

  this->buildM1(betaCoef);    // (σ dA/dt, A')
  this->buildCurl(alphaCoef); // (ν∇×A, ∇×A')
  this->buildGrad();
  b0 = new mfem::ParLinearForm(H1FESpace_);
  b1 = new mfem::ParLinearForm(HCurlFESpace_);
  A0 = new mfem::HypreParMatrix;
  A1 = new mfem::HypreParMatrix;
  A01 = new mfem::HypreParMatrix;
  A10 = new mfem::HypreParMatrix;
  blockA = new mfem::HypreParMatrix;

  X0 = new mfem::Vector;
  X1 = new mfem::Vector;
  B0 = new mfem::Vector;
  B1 = new mfem::Vector;

  mfem::Vector zero_vec(3);
  zero_vec = 0.0;
  mfem::VectorConstantCoefficient Zero_vec(zero_vec);
  mfem::ConstantCoefficient Zero(0.0);

  v_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  a_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[1]);

  v_.ProjectCoefficient(Zero);
  a_.ProjectCoefficient(Zero_vec);

  // aux
  e_.ProjectCoefficient(Zero_vec);
  b_.ProjectCoefficient(Zero_vec);
}

/*
This is the main computational code that computes dX/dt implicitly
where X is the state vector containing A and V

(M1+dt S1) E = WeakCurl^T B + Grad V
        S0 V = 0
*/
void AVSolver::ImplicitSolve(const double dt, const mfem::Vector &X,
                             mfem::Vector &dX_dt) {
  dX_dt = 0.0;
  dtCoef.constant = dt;
  v_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  a_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[1]);
  dv_.MakeRef(H1FESpace_, dX_dt, true_offsets[0]);
  da_.MakeRef(HCurlFESpace_, dX_dt, true_offsets[1]);

  _domain_properties.SetTime(this->GetTime());

  std::cout << "b0 start" << std::endl;
  mfem::ParGridFunction Phi_gf;
  Phi_gf.MakeRef(H1FESpace_, trueX->GetBlock(0), 0);

  mfem::Array<int> poisson_ess_tdof_list;
  Phi_gf = 1.0;
  *b0 = 0.0;
  std::cout << "b0 update start" << std::endl;
  // b0->Update(H1FESpace_, Phi_gf, 0);
  std::cout << "b0 update end" << std::endl;
  std::cout << "b0 bc start" << std::endl;

  _bc_map.applyIntegratedBCs(p_name, *b0, pmesh_);
  std::cout << "b0 ess start" << std::endl;
  _bc_map.applyEssentialBCs(p_name, poisson_ess_tdof_list, Phi_gf, pmesh_);
  std::cout << "b0 ess start" << std::endl;
  std::cout << "b0 ess start" << std::endl;
  std::cout << "b0 assemble" << std::endl;
  b0->Assemble();

  a0->EliminateEssentialBC(mfem::Array<int>({1, 1, 0}), Phi_gf, *b0);
  std::cout << "b0 ess startt" << std::endl;
  a0->Finalize();
  A0 = a0->ParallelAssemble();

  // _bc_map.applyEssentialBCs(p_name, poisson_ess_tdof_list, Phi_gf, pmesh_);
  std::cout << "b0 parallel assemble" << std::endl;

  // trueRhs->GetBlock(0).SyncAliasMemory(*trueRhs);

  // a0->EliminateEssentialBC(ess_bdr, Phi_gf, *b0);
  // a0->Finalize();

  // b0->ParallelAssemble(Phi_gf.GetTrueVector());
  std::cout << "b0 end" << std::endl;

  mfem::ParGridFunction J_gf(HCurlFESpace_);
  J_gf.MakeRef(HCurlFESpace_, trueX->GetBlock(1), 0);

  mfem::Array<int> ess_tdof_list;
  J_gf = 0.0;
  // b1->Update(HCurlFESpace_, J_gf, 0);
  *b1 = 0.0;
  curlCurl->MultTranspose(a_, *b1); // b1 = (ν∇×A, ∇×A')
  _bc_map.applyIntegratedBCs(u_name, *b1,
                             pmesh_); // b1 = (ν∇×A, ∇×A') +  <(ν∇×A) × n, A'>
  _bc_map.applyEssentialBCs(u_name, ess_tdof_list, J_gf, pmesh_);
  b1->Assemble();

  // b1->ParallelAssemble(J_gf.GetTrueVector());
  std::cout << "b1 end" << std::endl;

  // mVarf->EliminateEssentialBC(mark_bdr_attr_ess, trueX.GetBlock(0),
  // trueRhs.GetBlock(0)); bVarf->EliminateTrialDofs(mark_bdr_attr_ess,
  // trueX.GetBlock(0), trueRhs.GetBlock(1)); a01->Finalize();

  std::cout << "a1 start" << std::endl;

  mfem::Vector zero_vec(3);
  zero_vec = 0.0;

  if (a1 == NULL || fabs(dt - dt_A1) > 1.0e-12 * dt) {
    this->buildA1(betaCoef, dtAlphaCoef);
  }
  std::cout << "Assembled" << std::endl;
  std::cout << J_gf.Size() << std::endl;

  a1->EliminateEssentialBC(mfem::Array<int>({1, 1, 1}), J_gf, *b1);
  a1->Finalize();
  A1 = a1->ParallelAssemble();

  a01->EliminateTrialDofs(mfem::Array<int>({1, 1, 0}), J_gf, *b0);
  a10->EliminateTrialDofs(mfem::Array<int>({1, 1, 1}), Phi_gf, *b1);

  a01->Finalize();
  A01 = a01->ParallelAssemble();

  a10->Finalize();
  A10 = a10->ParallelAssemble();

  b0->ParallelAssemble(trueRhs->GetBlock(0));
  b1->ParallelAssemble(trueRhs->GetBlock(1));

  // A0->EliminateRowsCols(poisson_ess_tdof_list);
  // A1->EliminateRowsCols(ess_tdof_list);
  // A01->EliminateCols(ess_tdof_list);
  // A01->EliminateRows(poisson_ess_tdof_list);
  // A10->EliminateCols(poisson_ess_tdof_list);
  // A10->EliminateRows(ess_tdof_list);

  mfem::Array2D<mfem::HypreParMatrix *> hBlocks(2, 2);
  hBlocks = NULL;
  hBlocks(0, 0) = A0;
  hBlocks(0, 1) = A01;
  hBlocks(1, 0) = A10;
  hBlocks(1, 1) = A1;

  blockA = mfem::HypreParMatrixFromBlocks(hBlocks);
  blockA->Mult(*trueRhs, *trueX);

  // *trueX.Print();
  true_offsets.Print();
  block_trueOffsets.Print();
  std::cout << "trueX = " << trueX->Size() << std::endl;
  ;
  std::cout << "true RHS = " << trueRhs->Size() << std::endl;
  ;

  std::cout << "da = " << da_.Size() << std::endl;
  std::cout << "v = " << v_.Size() << std::endl;

  // solver.Mult(*trueRhs, *trueX);
  v_.Distribute(&(trueX->GetBlock(0)));
  da_.Distribute(&(trueX->GetBlock(1)));

  // A1 = HypreParMatrixFromBlocks(hBlocks);

  // if (ams_a1 == NULL) {
  //   mfem::ParFiniteElementSpace *prec_fespace =
  //       (a1->StaticCondensationIsEnabled() ? a1->SCParFESpace()
  //                                          : HCurlFESpace_);
  //   ams_a1 = new mfem::HypreAMS(*A1, prec_fespace);
  // }
  // if (pcg_a1 == NULL) {
  //   pcg_a1 = new mfem::HyprePCG(*A1);
  //   pcg_a1->SetTol(1.0e-9);
  //   pcg_a1->SetMaxIter(1000);
  //   pcg_a1->SetPrintLevel(0);
  //   pcg_a1->SetPreconditioner(*ams_a1);
  // }
  // // solve the system
  // // dE = (A1)^-1 [-S1 E]
  // pcg_a1->Mult(*B1, *X1);

  // a1->RecoverFEMSolution(*X1, *b1, e_);

  // ///

  // if (blockAV == NULL) {
  //   blockAV = new mfem::BlockOperator(block_trueOffsets);
  //   blockAV->SetBlock(0, 0, a0->ParallelAssemble());
  //   blockAV->SetBlock(0, 1, a01->ParallelAssemble());
  //   blockAV->SetBlock(1, 0, a10->ParallelAssemble());
  //   blockAV->SetBlock(1, 1, a1->ParallelAssemble());
  // }

  // if (amg_a0 == NULL) {
  //   amg_a0 = new mfem::HypreBoomerAMG(*A0);
  // }
  // if (ams_a1 == NULL) {
  //   mfem::ParFiniteElementSpace *prec_fespace =
  //       (a1->StaticCondensationIsEnabled() ? a1->SCParFESpace()
  //                                          : HCurlFESpace_);
  //   ams_a1 = new mfem::HypreAMS(*A1, prec_fespace);
  // }
  // // OperatorHandle ;
  // if (blockAVPr = NULL) {
  //   blockAVPr = new mfem::BlockDiagonalPreconditioner(block_trueOffsets);
  //   blockAVPr->SetDiagonalBlock(0, amg_a0);
  //   blockAVPr->SetDiagonalBlock(1, ams_a1);
  // }

  // if (pcg_a0 == NULL) {
  //   mfem::MINRESSolver solver(MPI_COMM_WORLD);
  //   solver.SetAbsTol(1.0e-9);
  //   //  solver.SetRelTol(rtol);
  //   solver.SetMaxIter(1000);
  //   solver.SetOperator(*blockAV);
  //   // solver.SetPreconditioner(*blockAVPr);
  //   //  solver.SetPrintLevel(verbose);
  //   solver.Mult(*trueRhs, *trueX);
  //   da_.Distribute(&(trueX->GetBlock(0)));
  //   v_.Distribute(&(trueX->GetBlock(1)));

  //   // pcg_a0 = new mfem::HyprePCG(MPI_COMM_WORLD);
  //   // pcg_a0->SetTol(1.0e-9);
  //   // pcg_a0->SetMaxIter(1000);
  //   // pcg_a0->SetPrintLevel(0);
  //   // pcg_a0->SetPreconditioner(dynamic_cast<mfem::HypreSolver
  //   *>(blockAVPr));
  //   // pcg_a0->SetOperator(*blockAV);
  // }

  // pcg_a0->Mult(*B0, *X0);
  // solver.Mult(trueRhs, trueX);

  // //  u->MakeTRef(R_space, x.GetBlock(0), 0);
  // //  p->MakeTRef(W_space, x.GetBlock(1), 0);
  // da_->Distribute(&(trueX.GetBlock(0)));
  // v_->Distribute(&(trueX.GetBlock(1)));

  // if (amg_a0 == NULL) {
  //   amg_a0 = new mfem::HypreBoomerAMG(*A0);
  // }
  // if (pcg_a0 == NULL) {
  //   pcg_a0 = new mfem::HyprePCG(*A0);
  //   pcg_a0->SetTol(1.0e-9);
  //   pcg_a0->SetMaxIter(1000);
  //   pcg_a0->SetPrintLevel(0);
  //   pcg_a0->SetPreconditioner(*amg_a0);
  //   blockAVPr->SetDiagonalBlock(0, pcg_a0);
  // }

  // We only need to create the solver and preconditioner once
  // if (ams_a1 == NULL) {
  //   mfem::ParFiniteElementSpace *prec_fespace =
  //       (a1->StaticCondensationIsEnabled() ? a1->SCParFESpace()
  //                                          : HCurlFESpace_);
  //   ams_a1 = new mfem::HypreAMS(*A1, prec_fespace);
  // }
  // if (pcg_a1 == NULL) {
  //   pcg_a1 = new mfem::HyprePCG(*A1);
  //   pcg_a1->SetTol(1.0e-9);
  //   pcg_a1->SetMaxIter(1000);
  //   pcg_a1->SetPrintLevel(0);
  //   pcg_a1->SetPreconditioner(*ams_a1);
  //   blockAVPr->SetDiagonalBlock(0, pcg_a1);
  // }

  // form the Laplacian and solve it

  // a0->FormLinearSystem(poisson_ess_tdof_list, Phi_gf, *b0, *A0, *X0, *B0);

  // if (amg_a0 == NULL) {
  //   amg_a0 = new mfem::HypreBoomerAMG(*A0);
  // }
  // if (pcg_a0 == NULL) {
  //   pcg_a0 = new mfem::HyprePCG(*A0);
  //   pcg_a0->SetTol(1.0e-9);
  //   pcg_a0->SetMaxIter(1000);
  //   pcg_a0->SetPrintLevel(0);
  //   pcg_a0->SetPreconditioner(*amg_a0);
  // }
  // pcg "Mult" operation is a solve
  // X0 = A0^-1 * B0

  // "undo" the static condensation saving result in grid function dP
  // a0->RecoverFEMSolution(*X0, *b0, v_);
  // dv_ = 0.0;
  // pcg_a1->Mult(*B1, *X1);

  // a1->RecoverFEMSolution(*X1, *b1, da_);

  //////////////////////////////////////////////////////////////////////////////
  // (σ(dA/dt + ∇ V), A') + (ν∇×A, ∇×A') - <(ν∇×A) × n, A'> - (J0, A') = 0

  // use a_ as a temporary, E = Grad P
  // b1 = -dt * Grad V

  // mfem::ParGridFunction J_gf(HCurlFESpace_);
  // mfem::Array<int> ess_tdof_list;
  // J_gf = 0.0;
  // _bc_map.applyEssentialBCs(u_name, ess_tdof_list, J_gf, pmesh_);
  // if (a1 == NULL || fabs(dt - dt_A1) > 1.0e-12 * dt) {
  //   this->buildA1(betaCoef, dtAlphaCoef);
  // }
  // // a1 = (σ dA/dt, A') + dt (ν∇×dA/dt, ∇×A')
  // // TODO: add  (σ ∇ V, A')
  // a1->FormLinearSystem(ess_tdof_list, J_gf, *b1, *A1, *X1, *B1);

  // // We only need to create the solver and preconditioner once
  // if (ams_a1 == NULL) {
  //   mfem::ParFiniteElementSpace *prec_fespace =
  //       (a1->StaticCondensationIsEnabled() ? a1->SCParFESpace()
  //                                          : HCurlFESpace_);
  //   ams_a1 = new mfem::HypreAMS(*A1, prec_fespace);
  // }
  // if (pcg_a1 == NULL) {
  //   pcg_a1 = new mfem::HyprePCG(*A1);
  //   pcg_a1->SetTol(1.0e-9);
  //   pcg_a1->SetMaxIter(1000);
  //   pcg_a1->SetPrintLevel(0);
  //   pcg_a1->SetPreconditioner(*ams_a1);
  // }
  // // solve the system
  // // dE = (A1)^-1 [-S1 E]
  // pcg_a1->Mult(*B1, *X1);

  // a1->RecoverFEMSolution(*X1, *b1, da_);

  // Auxiliary calculations
  // Compute B = Curl(A)
  curl->Mult(a_, b_);
  // Compute e = -(dA/dt + ∇ V);
  e_ = da_;
  grad->AddMult(v_, e_, 1.0);
  e_ *= -1.0;
}

void AVSolver::buildA1(mfem::Coefficient *Sigma, mfem::Coefficient *DtMuInv) {
  if (a1 != NULL) {
    delete a1;
  }
  // (σ(dA/dt, A') + dt (ν∇×dA/dt, ∇×A')
  // First create and assemble the bilinear form.  For now we assume the mesh
  // isn't moving, the materials are time independent, and dt is constant. So
  // we only need to do this once.

  a1 = new mfem::ParBilinearForm(HCurlFESpace_);
  a1->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*Sigma));
  a1->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*DtMuInv));
  a1->Assemble();

  // Don't finalize or parallel assemble this is done in FormLinearSystem.

  dt_A1 = dtCoef.constant;
}

void AVSolver::buildM1(mfem::Coefficient *Sigma) {
  if (m1 != NULL) {
    delete m1;
  }

  m1 = new mfem::ParBilinearForm(HCurlFESpace_);
  m1->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*Sigma));
  m1->Assemble();

  // Don't finalize or parallel assemble this is done in FormLinearSystem.
}

void AVSolver::buildGrad() {
  if (grad != NULL) {
    delete grad;
  }

  grad = new mfem::ParDiscreteLinearOperator(H1FESpace_, HCurlFESpace_);
  grad->AddDomainInterpolator(new mfem::GradientInterpolator());
  grad->Assemble();

  // no ParallelAssemble since this will be applied to GridFunctions
}

void AVSolver::buildCurl(mfem::Coefficient *MuInv) {
  if (curl != NULL) {
    delete curl;
  }
  if (weakCurl != NULL) {
    delete weakCurl;
  }

  curl = new mfem::ParDiscreteLinearOperator(HCurlFESpace_, HDivFESpace_);
  curl->AddDomainInterpolator(new mfem::CurlInterpolator);
  curl->Assemble();

  weakCurl = new mfem::ParMixedBilinearForm(HCurlFESpace_, HDivFESpace_);
  weakCurl->AddDomainIntegrator(new mfem::VectorFECurlIntegrator(*MuInv));
  weakCurl->Assemble();

  curlCurl = new mfem::ParBilinearForm(HCurlFESpace_);
  curlCurl->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*MuInv));
  curlCurl->Assemble();

  // no ParallelAssemble since this will be applied to GridFunctions
}

void AVSolver::SetVariableNames() {
  p_name = "electric_potential";
  p_display_name = "Electric Scalar Potential (V)";

  u_name = "magnetic_vector_potential";
  u_display_name = "Magnetic Vector Potential (A)";

  v_name = "magnetic_flux_density";
  v_display_name = "Magnetic Flux Density (B)";

  e_name = "electric_field";
  e_display_name = "Electric Field (E)";
}

void AVSolver::SetMaterialCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  alphaCoef = new mfem::TransformedCoefficient(
      &oneCoef, domain_properties.scalar_property_map["magnetic_permeability"],
      fracFunc);
  betaCoef = domain_properties.scalar_property_map["electrical_conductivity"];
}

void AVSolver::SetSourceCoefficient(
    hephaestus::DomainProperties &domain_properties) {
  if (domain_properties.vector_property_map.find("source") !=
      domain_properties.vector_property_map.end()) {
    sourceVecCoef = domain_properties.vector_property_map["source"];
  }
}

void AVSolver::buildSource() {
  // Replace with class to calculate div free source from input
  // VectorCoefficient
  src_gf = new mfem::ParGridFunction(HCurlFESpace_);
  div_free_src_gf = new mfem::ParGridFunction(HCurlFESpace_);
  _variables.Register("source", div_free_src_gf, false);
  // int irOrder = H1FESpace_->GetElementTransformation(0)->OrderW() +
  //               2 * H1FESpace_->GetOrder();
  int irOrder = H1FESpace_->GetElementTransformation(0)->OrderW() + 2 * 2;
  divFreeProj = new mfem::common::DivergenceFreeProjector(
      *H1FESpace_, *HCurlFESpace_, irOrder, NULL, NULL, NULL);
  hCurlMass = new mfem::ParBilinearForm(HCurlFESpace_);
  hCurlMass->AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
  hCurlMass->Assemble();
}

void AVSolver::RegisterOutputFields(mfem::DataCollection *dc_) {
  dc_->SetMesh(pmesh_);
  for (auto var = _variables.begin(); var != _variables.end(); ++var) {
    dc_->RegisterField(var->first, var->second);
  }
}

void AVSolver::WriteConsoleSummary(double t, int it) {
  // Write a summary of the timestep to console.
  if (myid_ == 0) {
    std::cout << std::fixed;
    std::cout << "step " << std::setw(6) << it << ",\tt = " << std::setw(6)
              << std::setprecision(3) << t << std::endl;
  }
}

void AVSolver::WriteOutputFields(mfem::DataCollection *dc_, int it) {
  if (dc_) {
    dc_->SetCycle(it);
    dc_->SetTime(t);
    dc_->Save();
  }
}

void AVSolver::InitializeGLVis() {
  if (myid_ == 0) {
    std::cout << "Opening GLVis sockets." << std::endl;
  }

  for (auto var = _variables.begin(); var != _variables.end(); ++var) {
    socks_[var->first] = new mfem::socketstream;
    socks_[var->first]->precision(8);
  }

  if (myid_ == 0) {
    std::cout << "GLVis sockets open." << std::endl;
  }
}

void AVSolver::DisplayToGLVis() {
  char vishost[] = "localhost";
  int visport = 19916;

  int Wx = 0, Wy = 0;                 // window position
  int Ww = 350, Wh = 350;             // window size
  int offx = Ww + 10, offy = Wh + 45; // window offsets

  for (auto var = _variables.begin(); var != _variables.end(); ++var) {
    mfem::common::VisualizeField(*socks_[var->first], vishost, visport,
                                 *(var->second), (var->first).c_str(), Wx, Wy,
                                 Ww, Wh);
    Wx += offx;
  }
}

} // namespace hephaestus
