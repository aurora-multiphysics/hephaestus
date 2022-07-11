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

AVSolver::AVSolver(mfem::ParMesh &pmesh, int order, hephaestus::BCMap &bc_map,
                   hephaestus::DomainProperties &domain_properties)
    : myid_(0), num_procs_(1), pmesh_(&pmesh), _bc_map(bc_map),
      _domain_properties(domain_properties),
      H1FESpace_(
          new mfem::common::H1_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HCurlFESpace_(
          new mfem::common::ND_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HDivFESpace_(
          new mfem::common::RT_ParFESpace(&pmesh, order, pmesh.Dimension())),
      true_offsets(4), a1(NULL), amg_a0(NULL), pcg_a0(NULL), ams_a1(NULL),
      pcg_a1(NULL), m1(NULL), grad(NULL), curl(NULL), weakCurl(NULL),
      v_(mfem::ParGridFunction(H1FESpace_)),
      a_(mfem::ParGridFunction(HCurlFESpace_)),
      b_(mfem::ParGridFunction(HDivFESpace_)),
      dv_(mfem::ParGridFunction(H1FESpace_)),
      da_(mfem::ParGridFunction(HCurlFESpace_)),
      db_(mfem::ParGridFunction(HDivFESpace_)) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);

  true_offsets[0] = 0;
  true_offsets[1] = H1FESpace_->GetVSize();
  true_offsets[2] = HCurlFESpace_->GetVSize();
  true_offsets[3] = HDivFESpace_->GetVSize();
  true_offsets.PartialSum();

  this->height = true_offsets[3];
  this->width = true_offsets[3];
}

void AVSolver::Init(mfem::Vector &X) {
  // Define material property coefficients
  dtCoef = mfem::ConstantCoefficient(1.0);
  oneCoef = mfem::ConstantCoefficient(1.0);
  SetMaterialCoefficients(_domain_properties);
  dtAlphaCoef = new mfem::TransformedCoefficient(&dtCoef, alphaCoef, prodFunc);

  SetVariableNames();

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

  this->buildM1(betaCoef);    // (σ dA/dt, A')
  this->buildCurl(alphaCoef); // (ν∇×A, ∇×A')
  this->buildGrad();
  b0 = new mfem::ParLinearForm(H1FESpace_);
  b1 = new mfem::ParLinearForm(HCurlFESpace_);
  A0 = new mfem::HypreParMatrix;
  A1 = new mfem::HypreParMatrix;
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
  b_.MakeRef(HDivFESpace_, const_cast<mfem::Vector &>(X), true_offsets[2]);

  v_.ProjectCoefficient(Zero);
  a_.ProjectCoefficient(Zero_vec);
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
  b_.MakeRef(HDivFESpace_, const_cast<mfem::Vector &>(X), true_offsets[2]);

  dv_.MakeRef(H1FESpace_, dX_dt, true_offsets[0]);
  da_.MakeRef(HCurlFESpace_, dX_dt, true_offsets[1]);
  db_.MakeRef(HDivFESpace_, dX_dt, true_offsets[2]);

  _domain_properties.SetTime(this->GetTime());

  // form the Laplacian and solve it
  mfem::ParGridFunction Phi_gf(H1FESpace_);
  mfem::Array<int> poisson_ess_tdof_list;
  Phi_gf = 0.0;
  *b0 = 0.0;
  _bc_map.applyEssentialBCs(p_name, poisson_ess_tdof_list, Phi_gf, pmesh_);
  _bc_map.applyIntegratedBCs(p_name, *b0, pmesh_);
  b0->Assemble();
  a0->FormLinearSystem(poisson_ess_tdof_list, Phi_gf, *b0, *A0, *X0, *B0);

  if (amg_a0 == NULL) {
    amg_a0 = new mfem::HypreBoomerAMG(*A0);
  }
  if (pcg_a0 == NULL) {
    pcg_a0 = new mfem::HyprePCG(*A0);
    pcg_a0->SetTol(1.0e-9);
    pcg_a0->SetMaxIter(1000);
    pcg_a0->SetPrintLevel(0);
    pcg_a0->SetPreconditioner(*amg_a0);
  }
  // pcg "Mult" operation is a solve
  // X0 = A0^-1 * B0
  pcg_a0->Mult(*B0, *X0);

  // "undo" the static condensation saving result in grid function dP
  a0->RecoverFEMSolution(*X0, *b0, v_);
  dv_ = 0.0;
  //////////////////////////////////////////////////////////////////////////////
  // (σ(dA/dt + ∇ V), A') + (ν∇×A, ∇×A') - <(ν∇×A) × n, A'> - (J0, A') = 0

  // use a_ as a temporary, E = Grad P
  // b1 = -dt * Grad V
  curlCurl->MultTranspose(a_, *b1); // b1 = (ν∇×A, ∇×A')
  grad->Mult(v_, da_);
  m1->AddMult(da_, *b1, 1.0); // b1 = (ν∇×A, ∇×A') + (J0, A')
  _bc_map.applyIntegratedBCs(
      u_name, *b1,
      pmesh_); // b1 = (ν∇×A, ∇×A') + (J0, A') + <(ν∇×A) × n, A'>

  mfem::ParGridFunction J_gf(HCurlFESpace_);
  mfem::Array<int> ess_tdof_list;
  J_gf = 0.0;
  _bc_map.applyEssentialBCs(u_name, ess_tdof_list, J_gf, pmesh_);
  if (a1 == NULL || fabs(dt - dt_A1) > 1.0e-12 * dt) {
    this->buildA1(betaCoef, dtAlphaCoef);
  }
  // a1 = (σ dA/dt, A') + dt (ν∇×dA/dt, ∇×A')
  // TODO: add  (σ ∇ V, A')
  a1->FormLinearSystem(ess_tdof_list, J_gf, *b1, *A1, *X1, *B1);

  // We only need to create the solver and preconditioner once
  if (ams_a1 == NULL) {
    mfem::ParFiniteElementSpace *prec_fespace =
        (a1->StaticCondensationIsEnabled() ? a1->SCParFESpace()
                                           : HCurlFESpace_);
    ams_a1 = new mfem::HypreAMS(*A1, prec_fespace);
  }
  if (pcg_a1 == NULL) {
    pcg_a1 = new mfem::HyprePCG(*A1);
    pcg_a1->SetTol(1.0e-9);
    pcg_a1->SetMaxIter(1000);
    pcg_a1->SetPrintLevel(0);
    pcg_a1->SetPreconditioner(*ams_a1);
  }
  // solve the system
  // dE = (A1)^-1 [-S1 E]
  pcg_a1->Mult(*B1, *X1);

  a1->RecoverFEMSolution(*X1, *b1, da_);

  // Compute B = Curl(A)
  curl->Mult(a_, b_);
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
}

void AVSolver::SetMaterialCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  alphaCoef = new mfem::TransformedCoefficient(
      &oneCoef, domain_properties.scalar_property_map["magnetic_permeability"],
      fracFunc);
  betaCoef = domain_properties.scalar_property_map["electrical_conductivity"];
}

void AVSolver::RegisterOutputFields(mfem::DataCollection *dc_) {
  dc_->RegisterField(u_name, &a_);
  dc_->RegisterField(v_name, &b_);
  dc_->RegisterField(p_name, &v_);
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

  socks_[u_name] = new mfem::socketstream;
  socks_[u_name]->precision(8);

  socks_[v_name] = new mfem::socketstream;
  socks_[v_name]->precision(8);

  socks_[p_name] = new mfem::socketstream;
  socks_[p_name]->precision(8);

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

  mfem::common::VisualizeField(*socks_[u_name], vishost, visport, a_,
                               u_display_name.c_str(), Wx, Wy, Ww, Wh);
  Wx += offx;

  mfem::common::VisualizeField(*socks_[v_name], vishost, visport, b_,
                               v_display_name.c_str(), Wx, Wy, Ww, Wh);
  Wx += offx;

  mfem::common::VisualizeField(*socks_[p_name], vishost, visport, v_,
                               p_display_name.c_str(), Wx, Wy, Ww, Wh);
  Wx += offx;
}

} // namespace hephaestus
