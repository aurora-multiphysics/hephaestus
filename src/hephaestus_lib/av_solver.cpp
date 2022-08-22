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
    : myid_(0), num_procs_(1), order_(order), pmesh_(&pmesh),
      _variables(variables), _bc_map(bc_map),
      _domain_properties(domain_properties),
      H1FESpace_(
          new mfem::common::H1_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HCurlFESpace_(
          new mfem::common::ND_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HDivFESpace_(
          new mfem::common::RT_ParFESpace(&pmesh, order, pmesh.Dimension())),
      a0(NULL), a1(NULL), a01(NULL), a10(NULL), amg_a0(NULL), pcg_a0(NULL),
      ams_a0(NULL), pcg_a1(NULL), m1(NULL), negBetaCoef(NULL), grad(NULL),
      curl(NULL), curlCurl(NULL), sourceVecCoef(NULL), src_gf(NULL),
      div_free_src_gf(NULL), hCurlMass(NULL), divFreeProj(NULL),
      p_(mfem::ParGridFunction(H1FESpace_)),
      u_(mfem::ParGridFunction(HCurlFESpace_)),
      dp_(mfem::ParGridFunction(H1FESpace_)),
      du_(mfem::ParGridFunction(HCurlFESpace_)),
      e_(mfem::ParGridFunction(HCurlFESpace_)),
      b_(mfem::ParGridFunction(HDivFESpace_)) {
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
  _variables.Register(e_name, &e_, false);
  _variables.Register(b_name, &b_, false);

  // Define material property coefficients
  dtCoef = mfem::ConstantCoefficient(1.0);
  oneCoef = mfem::ConstantCoefficient(1.0);
  SetMaterialCoefficients(_domain_properties);
  dtAlphaCoef = new mfem::TransformedCoefficient(&dtCoef, alphaCoef, prodFunc);
  SetSourceCoefficient(_domain_properties);
  if (sourceVecCoef) {
    buildSource();
  }

  this->buildCurl(alphaCoef); // (α∇×u_{n}, ∇×u')
  this->buildGrad();          // (s0_{n+1}, u')
  b0 = new mfem::ParLinearForm(HCurlFESpace_);
  b1 = new mfem::ParLinearForm(H1FESpace_);
  b01 = new mfem::ParLinearForm(HCurlFESpace_);
  b10 = new mfem::ParLinearForm(H1FESpace_);

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

  //////////////////////////////////////////////////////////////////////////////
  mfem::BlockVector x(true_offsets), rhs(true_offsets);
  mfem::BlockVector trueX(block_trueOffsets), trueRhs(block_trueOffsets);
  trueX = 0.0;
  *b0 = 0.0;

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
  curlCurl->MultTranspose(u_, *b0);
  *b0 *= -1.0;

  mfem::ParGridFunction J_gf(HCurlFESpace_);
  mfem::Array<int> ess_tdof_list;
  J_gf = 0.0;
  _bc_map.applyEssentialBCs(u_name, ess_tdof_list, J_gf, pmesh_);
  _bc_map.applyIntegratedBCs(u_name, *b0, pmesh_);
  b0->Assemble();

  if (src_gf) {
    src_gf->ProjectCoefficient(*sourceVecCoef);
    // Compute the discretely divergence-free portion of src_gf
    // divFreeProj->Mult(*src_gf, *div_free_src_gf);
    // Compute the dual of div_free_src_gf
    hCurlMass->AddMult(*src_gf, *b0);
  }

  mfem::ParGridFunction Phi_gf(H1FESpace_);
  mfem::Array<int> poisson_ess_tdof_list;
  Phi_gf = 0.0;
  *b1 = 0.0;
  _bc_map.applyEssentialBCs(p_name, poisson_ess_tdof_list, Phi_gf, pmesh_);
  _bc_map.applyIntegratedBCs(p_name, *b1, pmesh_);
  b1->Assemble();

  if (a0 == NULL || a10 == NULL || fabs(dt - dt_A0) > 1.0e-12 * dt) {
    this->buildA0(betaCoef, dtAlphaCoef);
  }

  if (a1 == NULL || a01 == NULL || fabs(dt - dt_A1) > 1.0e-12 * dt) {
    this->buildA1(betaCoef);
  }

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

  trueX.GetBlock(0).SyncAliasMemory(trueX);
  trueX.GetBlock(1).SyncAliasMemory(trueX);

  trueRhs.GetBlock(0).SyncAliasMemory(trueRhs);
  trueRhs.GetBlock(1).SyncAliasMemory(trueRhs);

  mfem::Array2D<mfem::HypreParMatrix *> hBlocks(2, 2);
  hBlocks = NULL;
  hBlocks(0, 0) = A0;
  hBlocks(0, 1) = A10;
  hBlocks(1, 0) = A01;
  hBlocks(1, 1) = A1;
  blockA = mfem::HypreParMatrixFromBlocks(hBlocks);

  mfem::MUMPSSolver mumps;
  mumps.SetPrintLevel(0);
  mumps.SetMatrixSymType(mfem::MUMPSSolver::MatType::UNSYMMETRIC);
  mumps.SetOperator(*blockA);
  mumps.Mult(trueRhs, trueX);

  // We only need to create the solver and preconditioner once
  // mfem::BlockOperator B(true_offsets, true_offsets);
  // B.SetBlock(0, 0, A0);
  // B.SetBlock(0, 1, A10);
  // B.SetBlock(1, 0, A01);
  // B.SetBlock(1, 1, A1);
  // if (ams_a0 == NULL) {
  //   ams_a0 = new mfem::HypreAMS(*A0, HCurlFESpace_);
  //   ams_a0->SetSingularProblem();
  // }
  // if (amg_a0 == NULL) {
  //   amg_a0 = new mfem::HypreBoomerAMG(*A1);
  // }
  // mfem::BlockDiagonalPreconditioner P(block_trueOffsets);
  // P.SetDiagonalBlock(0, ams_a0);

  // mfem::CGSolver pcg(MPI_COMM_WORLD);
  // pcg.SetOperator(B);
  // pcg.SetPreconditioner(P);
  // pcg.SetRelTol(1e-4);
  // pcg.SetMaxIter(1000);
  // pcg.SetPrintLevel(1);
  // pcg.Mult(trueRhs, trueX);

  trueX.GetBlock(0).SyncAliasMemory(trueX);
  trueX.GetBlock(1).SyncAliasMemory(trueX);

  du_.Distribute(&(trueX.GetBlock(0)));
  p_.Distribute(&(trueX.GetBlock(1)));

  grad->Mult(p_, e_);
  e_ += du_;
  e_ *= -1.0;

  curl->Mult(u_, b_);
  curl->AddMult(du_, b_, dt);
}

void AVSolver::buildA0(mfem::Coefficient *betaCoef,
                       mfem::Coefficient *dtAlphaCoef) {
  if (a0 != NULL) {
    delete a0;
  }
  if (a10 != NULL) {
    delete a10;
  }

  // a0(dA/dt, dA'/dt) = (βdA/dt_{n+1}, dA'/dt) + (αdt∇×dA/dt_{n+1}, ∇×dA'/dt)
  a0 = new mfem::ParBilinearForm(HCurlFESpace_);
  a0->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*betaCoef));
  a0->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*dtAlphaCoef));
  a0->Assemble();

  // a10(V, dA'/dt) = (σ ∇ V, dA'/dt)
  a10 = new mfem::ParMixedBilinearForm(H1FESpace_, HCurlFESpace_);
  a10->AddDomainIntegrator(new mfem::MixedVectorGradientIntegrator(*betaCoef));
  a10->Assemble();

  dt_A0 = dtCoef.constant;
}

void AVSolver::buildA1(mfem::Coefficient *betaCoef) {
  if (a1 != NULL) {
    delete a1;
  }
  if (a01 != NULL) {
    delete a01;
  }
  if (negBetaCoef != NULL) {
    delete negBetaCoef;
  }

  negCoef = mfem::ConstantCoefficient(-1.0);
  negBetaCoef = new mfem::TransformedCoefficient(&negCoef, betaCoef, prodFunc);

  // a1(V, V') = (σ ∇ V, ∇ V')
  a1 = new mfem::ParBilinearForm(H1FESpace_);
  a1->AddDomainIntegrator(new mfem::DiffusionIntegrator(*negBetaCoef));
  a1->Assemble();

  // (σdA/dt, ∇ V')
  a01 = new mfem::ParMixedBilinearForm(HCurlFESpace_, H1FESpace_);
  a01->AddDomainIntegrator(
      new mfem::VectorFEWeakDivergenceIntegrator(*negBetaCoef));
  a01->Assemble();

  dt_A1 = dtCoef.constant;
}

void AVSolver::buildGrad() {
  // Discrete Grad operator
  if (grad != NULL) {
    delete grad;
  }
  grad = new mfem::ParDiscreteLinearOperator(H1FESpace_, HCurlFESpace_);
  grad->AddDomainInterpolator(new mfem::GradientInterpolator());
  grad->Assemble();
}

void AVSolver::buildCurl(mfem::Coefficient *MuInv) {
  if (curlCurl != NULL) {
    delete curlCurl;
  }
  curlCurl = new mfem::ParBilinearForm(HCurlFESpace_);
  curlCurl->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*MuInv));
  curlCurl->Assemble();

  // Discrete Curl operator
  if (curl != NULL) {
    delete curl;
  }
  curl = new mfem::ParDiscreteLinearOperator(HCurlFESpace_, HDivFESpace_);
  curl->AddDomainInterpolator(new mfem::CurlInterpolator());
  curl->Assemble();
}

void AVSolver::SetVariableNames() {
  p_name = "electric_potential";
  p_display_name = "Electric Scalar Potential";

  u_name = "magnetic_vector_potential";
  u_display_name = "Magnetic Vector Potential";

  e_name = "electric_field";
  e_display_name = "Electric Field";

  b_name = "magnetic_flux_density";
  b_display_name = "Magnetic Flux Density";
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
  _variables.Register("source", src_gf, false);
  int irOrder = H1FESpace_->GetElementTransformation(0)->OrderW() + 2 * order_;
  // divFreeProj = new mfem::common::DivergenceFreeProjector(
  //     *H1FESpace_, *HCurlFESpace_, irOrder, NULL, NULL, NULL);
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
