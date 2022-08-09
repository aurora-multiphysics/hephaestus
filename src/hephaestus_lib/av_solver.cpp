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
                   mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
                   hephaestus::BCMap &bc_map,
                   hephaestus::DomainProperties &domain_properties)
    : myid_(0), num_procs_(1), pmesh_(&pmesh), _variables(variables),
      _bc_map(bc_map), _domain_properties(domain_properties),
      H1FESpace_(
          new mfem::common::H1_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HCurlFESpace_(
          new mfem::common::ND_ParFESpace(&pmesh, order, pmesh.Dimension())),
      a0(NULL), a1(NULL), a01(NULL), a10(NULL), amg_a0(NULL), pcg_a0(NULL),
      ams_a0(NULL), pcg_a1(NULL), m1(NULL), grad(NULL), curl(NULL),
      curlCurl(NULL), sourceVecCoef(NULL), src_gf(NULL), div_free_src_gf(NULL),
      hCurlMass(NULL), divFreeProj(NULL), p_(mfem::ParGridFunction(H1FESpace_)),
      u_(mfem::ParGridFunction(HCurlFESpace_)),
      dp_(mfem::ParGridFunction(H1FESpace_)),
      du_(mfem::ParGridFunction(HCurlFESpace_)) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);

  true_offsets.SetSize(3);
  true_offsets[0] = 0;
  true_offsets[1] = HCurlFESpace_->GetVSize();
  true_offsets[2] = H1FESpace_->GetVSize();
  true_offsets.PartialSum();

  block_trueOffsets.SetSize(3);
  block_trueOffsets[0] = 0;
  block_trueOffsets[1] = HCurlFESpace_->TrueVSize();
  block_trueOffsets[2] = H1FESpace_->TrueVSize();
  block_trueOffsets.PartialSum();

  this->height = true_offsets[2];
  this->width = true_offsets[2];
}

void AVSolver::Init(mfem::Vector &X) {
  SetVariableNames();
  _variables.Register(u_name, &u_, false);
  _variables.Register(p_name, &p_, false);

  // Define material property coefficients
  dtCoef = mfem::ConstantCoefficient(1.0);
  oneCoef = mfem::ConstantCoefficient(1.0);
  negCoef = mfem::ConstantCoefficient(-1.0);
  SetMaterialCoefficients(_domain_properties);
  dtAlphaCoef = new mfem::TransformedCoefficient(&dtCoef, alphaCoef, prodFunc);
  negBetaCoef = new mfem::TransformedCoefficient(&negCoef, betaCoef, prodFunc);
  SetSourceCoefficient(_domain_properties);
  if (sourceVecCoef) {
    buildSource();
  }

  this->buildCurl(alphaCoef); // (α∇×u_{n}, ∇×u')
  this->buildGrad();          // (s0_{n+1}, u')
  b0 = new mfem::ParLinearForm(HCurlFESpace_);
  b1 = new mfem::ParLinearForm(H1FESpace_);
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

  u_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  p_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[1]);

  u_.ProjectCoefficient(Zero_vec);
  p_.ProjectCoefficient(Zero);
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
  dtCoef.constant = dt;

  u_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  p_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[1]);

  du_.MakeRef(HCurlFESpace_, dX_dt, true_offsets[0]);
  dp_.MakeRef(H1FESpace_, dX_dt, true_offsets[1]);

  _domain_properties.SetTime(this->GetTime());

  // // -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
  // // a0(p_{n+1}, p') = b0(p')
  // // a0(p, p') = (β ∇ p, ∇ p')
  // // b0(p') = <n.s0, p'>
  // mfem::ParGridFunction Phi_gf(H1FESpace_);
  // mfem::Array<int> poisson_ess_tdof_list;
  // Phi_gf = 0.0;
  // *b0 = 0.0;
  // _bc_map.applyEssentialBCs(p_name, poisson_ess_tdof_list, Phi_gf, pmesh_);
  // _bc_map.applyIntegratedBCs(p_name, *b0, pmesh_);
  // b0->Assemble();
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
  // // pcg "Mult" operation is a solve
  // // X0 = A0^-1 * B0
  // pcg_a0->Mult(*B0, *X0);

  // // "undo" the static condensation saving result in grid function dP
  // a0->RecoverFEMSolution(*X0, *b0, p_);
  dp_ = 0.0;
  //////////////////////////////////////////////////////////////////////////////
  // (α∇×u_{n}, ∇×u') + (αdt∇×du/dt_{n+1}, ∇×u') + (βdu/dt_{n+1}, u')
  // - (s0_{n+1}, u') - <(α∇×u_{n+1}) × n, u'> = 0
  mfem::BlockVector x(true_offsets), rhs(true_offsets);
  mfem::BlockVector trueX(block_trueOffsets), trueRhs(block_trueOffsets);
  trueX = 0.0;
  *b0 = 0.0;
  // b0->Update(HCurlFESpace_, rhs.GetBlock(0), 0);
  // a1(du/dt_{n+1}, u') = b1(u')
  // a1(u, u') = (βu, u') + (αdt∇×u, ∇×u')
  // b1(u') = (s0_{n+1}, u') - (α∇×u_{n}, ∇×u') + <(α∇×u_{n+1}) × n, u'>

  // (α∇×u_{n}, ∇×u')
  // v_ is a grid function but curlCurl is not parallel assembled so is OK
  curlCurl->MultTranspose(u_, *b0);
  *b0 *= -1.0;

  // use du_ as a temporary
  // (s0_{n+1}, u')
  // grad->Mult(p_, du_);
  // m1->AddMult(du_, *b1, 1.0);

  mfem::ParGridFunction J_gf(HCurlFESpace_);
  mfem::Array<int> ess_tdof_list;
  J_gf = 0.0;
  _bc_map.applyEssentialBCs(u_name, ess_tdof_list, J_gf, pmesh_);
  _bc_map.applyIntegratedBCs(u_name, *b0, pmesh_);
  b0->Assemble();

  if (src_gf) {
    src_gf->ProjectCoefficient(*sourceVecCoef);
    // Compute the discretely divergence-free portion of src_gf
    divFreeProj->Mult(*src_gf, *div_free_src_gf);
    // Compute the dual of div_free_src_gf
    hCurlMass->AddMult(*div_free_src_gf, *b0);
  }

  // b0->SyncAliasMemory(rhs);
  // trueRhs.GetBlock(0).SyncAliasMemory(trueRhs);

  // a1(du/dt_{n+1}, u') = (βdu/dt_{n+1}, u') + (αdt∇×du/dt_{n+1}, ∇×u')
  // if (a0 == NULL || fabs(dt - dt_A0) > 1.0e-12 * dt) {
  this->buildA0(betaCoef, dtAlphaCoef);
  // }
  // a0->EliminateEssentialBC(mfem::Array<int>({1, 1, 1}), J_gf, *b0);

  // a0->EliminateEssentialBC(mfem::Array<int>({1, 1, 1}), J_gf, *b0);
  // a0->Finalize();
  // A0 = a0->ParallelAssemble();

  // a0->FormSystemMatrix(ess_tdof_list, *A0);
  // A0->EliminateRowsCols(ess_tdof_list);

  // a0->FormLinearSystem(ess_tdof_list, J_gf, *b0, *A0, *X0, *B0);

  mfem::ParGridFunction Phi_gf(H1FESpace_);
  mfem::Array<int> poisson_ess_tdof_list;
  Phi_gf = 0.0;
  *b1 = 0.0;
  _bc_map.applyEssentialBCs(p_name, poisson_ess_tdof_list, Phi_gf, pmesh_);
  _bc_map.applyIntegratedBCs(p_name, *b1, pmesh_);
  b1->Assemble();

  // a1->EliminateEssentialBC(mfem::Array<int>({1, 1, 1}), Phi_gf, *b1);
  // a1->Finalize();
  // A1 = a1->ParallelAssemble();

  delete a10;
  delete a01;
  delete a1;

  // a1(V, V') = (σ ∇ V, ∇ V')
  a1 = new mfem::ParBilinearForm(H1FESpace_);
  a1->AddDomainIntegrator(new mfem::DiffusionIntegrator(*negBetaCoef));
  a1->Assemble();

  // (σdA/dt , ∇ V')
  a01 = new mfem::ParMixedBilinearForm(HCurlFESpace_, H1FESpace_);
  a01->AddDomainIntegrator(
      new mfem::VectorFEWeakDivergenceIntegrator(*negBetaCoef));
  a01->Assemble();

  // (σ ∇ V, dA'/dt)
  a10 = new mfem::ParMixedBilinearForm(H1FESpace_, HCurlFESpace_);
  a10->AddDomainIntegrator(new mfem::MixedVectorGradientIntegrator(*betaCoef));
  a10->Assemble();

  // a01->EliminateTrialDofs(mfem::Array<int>({1, 1, 1}), J_gf, *b1);
  // a01->Finalize();
  // A01 = a01->ParallelAssemble();

  // a10->EliminateTrialDofs(mfem::Array<int>({1, 1, 1}), Phi_gf, *b0);
  // a10->Finalize();
  // A10 = a10->ParallelAssemble();

  mfem::ParLinearForm *b01 = new mfem::ParLinearForm(HCurlFESpace_);
  mfem::ParLinearForm *b10 = new mfem::ParLinearForm(H1FESpace_);

  *b01 = 0.0;
  *b10 = 0.0;

  a0->FormLinearSystem(ess_tdof_list, J_gf, *b0, *A0, *X0, *B0);
  trueX.GetBlock(0) = *X0;
  trueRhs.GetBlock(0) = *B0;
  a1->FormLinearSystem(poisson_ess_tdof_list, Phi_gf, *b1, *A1, *X1, *B1);
  trueX.GetBlock(1) = *X1;
  trueRhs.GetBlock(1) = *B1;

  a01->FormRectangularLinearSystem(ess_tdof_list, poisson_ess_tdof_list, J_gf,
                                   *b10, *A01, *X0, *B1);

  a10->FormRectangularLinearSystem(poisson_ess_tdof_list, ess_tdof_list, Phi_gf,
                                   *b01, *A10, *X1, *B0);
  trueRhs.GetBlock(1) += *B1;

  trueRhs.GetBlock(0) += *B0;

  // b0->ParallelAssemble(trueRhs.GetBlock(0));
  trueX.GetBlock(0).SyncAliasMemory(trueX);
  trueX.GetBlock(1).SyncAliasMemory(trueX);

  trueRhs.GetBlock(0).SyncAliasMemory(trueRhs);
  trueRhs.GetBlock(1).SyncAliasMemory(trueRhs);

  // A0->EliminateRowsCols(ess_tdof_list, trueX.GetBlock(0),
  // trueRhs.GetBlock(0)); A1->EliminateRowsCols(poisson_ess_tdof_list,
  // trueX.GetBlock(1),
  //                       trueRhs.GetBlock(1));
  // A01->EliminateCols(ess_tdof_list);
  // A01->EliminateRows(poisson_ess_tdof_list);
  // A10->EliminateRows(ess_tdof_list);
  // A10->EliminateCols(poisson_ess_tdof_list);

  mfem::Array2D<mfem::HypreParMatrix *> hBlocks(2, 2);
  hBlocks = NULL;
  hBlocks(0, 0) = A0;
  hBlocks(0, 1) = A10;
  hBlocks(1, 0) = A01;
  hBlocks(1, 1) = A1;

  blockA = mfem::HypreParMatrixFromBlocks(hBlocks);

  // We only need to create the solver and preconditioner once
  // if (ams_a0 == NULL) {
  //   ams_a0 = new mfem::HypreAMS(*blockA, HCurlFESpace_);
  //   ams_a0->SetSingularProblem();
  // }
  mfem::MUMPSSolver mumps;
  mumps.SetPrintLevel(0);
  mumps.SetMatrixSymType(mfem::MUMPSSolver::MatType::UNSYMMETRIC);
  mumps.SetOperator(*blockA);
  mumps.Mult(trueRhs, trueX);

  // if (pcg_a0 == NULL) {
  //   pcg_a0 = new mfem::HyprePCG(*blockA);
  //   pcg_a0->SetTol(1.0e-16);
  //   pcg_a0->SetMaxIter(1000);
  //   pcg_a0->SetPrintLevel(0);
  //   // pcg_a0->SetPreconditioner(*ams_a0);
  // }
  // // solve the system
  // pcg_a0->Mult(trueRhs, trueX);

  trueX.GetBlock(0).SyncAliasMemory(trueX);
  trueX.GetBlock(1).SyncAliasMemory(trueX);

  du_.Distribute(&(trueX.GetBlock(0)));
  // du_ = trueX.GetBlock(0);
  //  p->MakeRef(W_space, x.GetBlock(1), 0);
  // du_.Distribute(&(trueX.GetBlock(0)));
  p_.Distribute(&(trueX.GetBlock(1)));

  // du_ = trueX;
  // a0->RecoverFEMSolution((trueX.GetBlock(0)), *b0, du_);
}

void AVSolver::buildA0(mfem::Coefficient *Sigma, mfem::Coefficient *DtMuInv) {
  if (a0 != NULL) {
    delete a0;
  }

  // First create and assemble the bilinear form.  For now we assume the mesh
  // isn't moving, the materials are time independent, and dt is constant. So
  // we only need to do this once.

  a0 = new mfem::ParBilinearForm(HCurlFESpace_);
  a0->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*Sigma));
  a0->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*DtMuInv));
  a0->Assemble();

  // Don't finalize or parallel assemble this is done in FormLinearSystem.

  dt_A0 = dtCoef.constant;
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
  if (curlCurl != NULL) {
    delete curlCurl;
  }

  curlCurl = new mfem::ParBilinearForm(HCurlFESpace_);
  curlCurl->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*MuInv));
  curlCurl->Assemble();

  // no ParallelAssemble since this will be applied to GridFunctions
}

void AVSolver::SetVariableNames() {
  p_name = "electric_potential";
  p_display_name = "Electric Scalar Potential";

  u_name = "magnetic_vector_potential";
  u_display_name = "Magnetic Vector Potential";
}

void AVSolver::SetMaterialCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  // if (domain_properties.scalar_property_map.count("magnetic_permeability") ==
  //     0) {
  //   domain_properties.scalar_property_map["magnetic_permeability"] =
  //       new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
  //           std::string("magnetic_permeability")));
  // }
  // if (domain_properties.scalar_property_map.count("electrical_conductivity")
  // ==
  //     0) {
  //   domain_properties.scalar_property_map["electrical_conductivity"] =
  //       new mfem::PWCoefficient(domain_properties.getGlobalScalarProperty(
  //           std::string("electrical_conductivity")));
  // }
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
  // _variables.Register("source", src_gf, false);
  // _variables.Register("sourcedivfree", div_free_src_gf, false);
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
