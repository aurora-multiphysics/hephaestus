// Solves the equations
// ∇⋅s0 = 0
// ∇×(α∇×u) - βu = s0
// dv/dt = -∇×u

// where
// s0 ∈ H(div) source field
// v ∈ H(div)
// u ∈ H(curl)
// p ∈ H1

// Weak form (Space discretisation)
// -(s0, ∇ p') + <n.s0, p'> = 0
// (α∇×u, ∇×u') - (βu, u') - (s0, u') - <(α∇×u) × n, u'> = 0
// (dv/dt, v') + (∇×u, v') = 0

// Time discretisation using implicit scheme:
// Unknowns
// s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
// dv/dt_{n+1} ∈ H(div)
// u_{n+1} ∈ H(curl)
// p_{n+1} ∈ H1

// Fully discretised equations
// -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
// (αp_{n}, ∇×u') - (αdt∇×u_{n+1}, ∇×u') - (βu_{n+1}, u') - (s0_{n+1}, u') -
// <(α∇×u_{n+1}) × n, u'> = 0
// (dv/dt_{n+1}, v') + (∇×u_{n+1}, v') = 0
// using
// p_{n+1} = p_{n} + dt dv/dt_{n+1} = p_{n} - dt ∇×u_{n+1}

// Rewritten as
// a0(p_{n+1}, p') = b0(p')
// a1(u_{n+1}, u') = b1(u')
// dv/dt_{n+1} = -∇×u

// where
// a0(p, p') = (β ∇ p, ∇ p')
// b0(p') = <n.s0, p'>
// a1(u, u') = (βu, u') + (αdt∇×u, ∇×u')
// b1(u') = (s0_{n+1}, u') + (αp_{n}, ∇×u') + <(αdt∇×u_{n+1}) × n, u'>

#include "hcurl_solver.hpp"

namespace hephaestus {
double prodFunc(double a, double b) { return a * b; }
double fracFunc(double a, double b) { return a / b; }

HCurlSolver::HCurlSolver(mfem::ParMesh &pmesh, int order,
                         hephaestus::BCMap &bc_map,
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
      p_(mfem::ParGridFunction(H1FESpace_)),
      u_(mfem::ParGridFunction(HCurlFESpace_)),
      v_(mfem::ParGridFunction(HDivFESpace_)),
      dp_(mfem::ParGridFunction(H1FESpace_)),
      du_(mfem::ParGridFunction(HCurlFESpace_)),
      dv_(mfem::ParGridFunction(HDivFESpace_)) {
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

void HCurlSolver::Init(mfem::Vector &X) {
  // Define material property coefficients
  dtCoef = mfem::ConstantCoefficient(1.0);
  oneCoef = mfem::ConstantCoefficient(1.0);
  SetMaterialCoefficients(_domain_properties);
  dtAlphaCoef = new mfem::TransformedCoefficient(&dtCoef, alphaCoef, prodFunc);

  SetVariableNames();

  // a0(p, p') = (β ∇ p, ∇ p')
  a0 = new mfem::ParBilinearForm(H1FESpace_);
  a0->AddDomainIntegrator(new mfem::DiffusionIntegrator(*betaCoef));
  a0->Assemble();

  this->buildM1(betaCoef);    // (βu, u')
  this->buildCurl(alphaCoef); // (αp_{n}, ∇×u')
  this->buildGrad();          // (s0_{n+1}, u')
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

  p_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  u_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[1]);
  v_.MakeRef(HDivFESpace_, const_cast<mfem::Vector &>(X), true_offsets[2]);

  p_.ProjectCoefficient(Zero);
  u_.ProjectCoefficient(Zero_vec);
  v_.ProjectCoefficient(Zero_vec);
}

/*
This is the main computational code that computes dX/dt implicitly
where X is the state vector containing p, u and v.

Unknowns
s0_{n+1} ∈ H(div) source field, where s0 = -β∇p
dv/dt_{n+1} ∈ H(div)
u_{n+1} ∈ H(curl)
p_{n+1} ∈ H1

Fully discretised equations
-(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
(αp_{n}, ∇×u') - (αdt∇×u_{n+1}, ∇×u') - (βu_{n+1}, u') - (s0_{n+1}, u') -
<(α∇×u_{n+1}) × n, u'> = 0
(dv/dt_{n+1}, v') + (∇×u_{n+1}, v') = 0
using
p_{n+1} = p_{n} + dt dv/dt_{n+1} = p_{n} - dt ∇×u_{n+1}
*/
void HCurlSolver::ImplicitSolve(const double dt, const mfem::Vector &X,
                                mfem::Vector &dX_dt) {
  dX_dt = 0.0;
  dtCoef.constant = dt;

  p_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  u_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[1]);
  v_.MakeRef(HDivFESpace_, const_cast<mfem::Vector &>(X), true_offsets[2]);

  dp_.MakeRef(H1FESpace_, dX_dt, true_offsets[0]);
  du_.MakeRef(HCurlFESpace_, dX_dt, true_offsets[1]);
  dv_.MakeRef(HDivFESpace_, dX_dt, true_offsets[2]);

  _domain_properties.SetTime(this->GetTime());

  // -(s0_{n+1}, ∇ p') + <n.s0_{n+1}, p'> = 0
  // a0(p_{n+1}, p') = b0(p')
  // a0(p, p') = (β ∇ p, ∇ p')
  // b0(p') = <n.s0, p'>
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
  a0->RecoverFEMSolution(*X0, *b0, p_);
  dp_ = 0.0;
  //////////////////////////////////////////////////////////////////////////////
  // (αp_{n}, ∇×u') - (αdt∇×u_{n+1}, ∇×u') - (βu_{n+1}, u') - (s0_{n+1}, u') -
  // <(α∇×u_{n+1}) × n, u'> = 0

  // a1(u_{n+1}, u') = b1(u')
  // a1(u, u') = (βu, u') + (αdt∇×u, ∇×u')
  // b1(u') = (s0_{n+1}, u') + (αp_{n}, ∇×u') + <(αdt∇×u_{n+1}) × n, u'>

  // (αp_{n}, ∇×u')
  // v_ is a grid function but weakCurl is not parallel assembled so is OK
  weakCurl->MultTranspose(v_, *b1);

  // use u_ as a temporary, E = Grad v
  // (s0_{n+1}, u')
  grad->Mult(p_, u_);
  m1->AddMult(u_, *b1, 1.0);

  mfem::ParGridFunction J_gf(HCurlFESpace_);
  mfem::Array<int> ess_tdof_list;
  J_gf = 0.0;
  _bc_map.applyEssentialBCs(u_name, ess_tdof_list, J_gf, pmesh_);
  _bc_map.applyIntegratedBCs(u_name, *b1, pmesh_);
  if (a1 == NULL || fabs(dt - dt_A1) > 1.0e-12 * dt) {
    this->buildA1(betaCoef, dtAlphaCoef);
  }
  a1->FormLinearSystem(ess_tdof_list, J_gf, *b1, *A1, *X1, *B1);

  // We only need to create the solver and preconditioner once
  if (ams_a1 == NULL) {
    ams_a1 = new mfem::HypreAMS(*A1, HCurlFESpace_);
  }
  if (pcg_a1 == NULL) {
    pcg_a1 = new mfem::HyprePCG(*A1);
    pcg_a1->SetTol(1.0e-9);
    pcg_a1->SetMaxIter(1000);
    pcg_a1->SetPrintLevel(0);
    pcg_a1->SetPreconditioner(*ams_a1);
  }
  // solve the system
  pcg_a1->Mult(*B1, *X1);

  a1->RecoverFEMSolution(*X1, *b1, u_);
  du_ = 0.0;

  // Subtract off contribution from source
  grad->AddMult(p_, u_, -1.0);

  // dv/dt_{n+1} = -∇×u
  // note curl maps GF to GF
  curl->Mult(u_, dv_);
  dv_ *= -1.0;
}

void HCurlSolver::buildA1(mfem::Coefficient *Sigma,
                          mfem::Coefficient *DtMuInv) {
  if (a1 != NULL) {
    delete a1;
  }

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

void HCurlSolver::buildM1(mfem::Coefficient *Sigma) {
  if (m1 != NULL) {
    delete m1;
  }

  m1 = new mfem::ParBilinearForm(HCurlFESpace_);
  m1->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*Sigma));
  m1->Assemble();

  // Don't finalize or parallel assemble this is done in FormLinearSystem.
}

void HCurlSolver::buildGrad() {
  if (grad != NULL) {
    delete grad;
  }

  grad = new mfem::ParDiscreteLinearOperator(H1FESpace_, HCurlFESpace_);
  grad->AddDomainInterpolator(new mfem::GradientInterpolator());
  grad->Assemble();

  // no ParallelAssemble since this will be applied to GridFunctions
}

void HCurlSolver::buildCurl(mfem::Coefficient *MuInv) {
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

  // no ParallelAssemble since this will be applied to GridFunctions
}

void HCurlSolver::SetVariableNames() {
  p_name = "scalar_potential";
  p_display_name = "Scalar Potential";

  u_name = "h_curl_var";
  u_display_name = "H(Curl) variable";

  v_name = "h_div_var";
  v_display_name = "H(Div) variable";
}

void HCurlSolver::SetMaterialCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  if (domain_properties.scalar_property_map.count("alpha") == 0) {
    domain_properties.scalar_property_map["alpha"] = new mfem::PWCoefficient(
        domain_properties.getGlobalScalarProperty(std::string("alpha")));
  }
  if (domain_properties.scalar_property_map.count("beta") == 0) {
    domain_properties.scalar_property_map["beta"] = new mfem::PWCoefficient(
        domain_properties.getGlobalScalarProperty(std::string("beta")));
  }
  alphaCoef = domain_properties.scalar_property_map["alpha"];
  betaCoef = domain_properties.scalar_property_map["beta"];
}

void HCurlSolver::RegisterOutputFields(mfem::DataCollection *dc_) {
  dc_->RegisterField(u_name, &u_);
  dc_->RegisterField(v_name, &v_);
  dc_->RegisterField(p_name, &p_);
}

void HCurlSolver::WriteConsoleSummary(double t, int it) {
  // Write a summary of the timestep to console.
  if (myid_ == 0) {
    std::cout << std::fixed;
    std::cout << "step " << std::setw(6) << it << ",\tt = " << std::setw(6)
              << std::setprecision(3) << t << std::endl;
  }
}

void HCurlSolver::WriteOutputFields(mfem::DataCollection *dc_, int it) {
  if (dc_) {
    dc_->SetCycle(it);
    dc_->SetTime(t);
    dc_->Save();
  }
}

void HCurlSolver::InitializeGLVis() {
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

void HCurlSolver::DisplayToGLVis() {
  char vishost[] = "localhost";
  int visport = 19916;

  int Wx = 0, Wy = 0;                 // window position
  int Ww = 350, Wh = 350;             // window size
  int offx = Ww + 10, offy = Wh + 45; // window offsets

  mfem::common::VisualizeField(*socks_[u_name], vishost, visport, u_,
                               u_display_name.c_str(), Wx, Wy, Ww, Wh);
  Wx += offx;

  mfem::common::VisualizeField(*socks_[v_name], vishost, visport, v_,
                               v_display_name.c_str(), Wx, Wy, Ww, Wh);
  Wx += offx;

  mfem::common::VisualizeField(*socks_[p_name], vishost, visport, p_,
                               p_display_name.c_str(), Wx, Wy, Ww, Wh);
  Wx += offx;
}

} // namespace hephaestus
