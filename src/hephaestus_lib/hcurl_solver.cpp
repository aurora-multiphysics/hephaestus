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

#include "hcurl_solver.hpp"

namespace hephaestus {

HCurlSolver::HCurlSolver(mfem::ParMesh &pmesh, int order,
                         hephaestus::BCMap &bc_map,
                         hephaestus::DomainProperties &domain_properties)
    : myid_(0), num_procs_(1), pmesh_(&pmesh), _bc_map(bc_map),
      _domain_properties(domain_properties),
      H1FESpace_(
          new mfem::common::H1_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HCurlFESpace_(
          new mfem::common::ND_ParFESpace(&pmesh, order, pmesh.Dimension())),
      a1(NULL), amg_a0(NULL), pcg_a0(NULL), ams_a1(NULL), pcg_a1(NULL),
      m1(NULL), grad(NULL), curl(NULL), curlCurl(NULL), sourceVecCoef(NULL),
      src_gf(NULL), div_free_src_gf(NULL), hCurlMass(NULL), divFreeProj(NULL),
      p_(mfem::ParGridFunction(H1FESpace_)),
      u_(mfem::ParGridFunction(HCurlFESpace_)),
      dp_(mfem::ParGridFunction(H1FESpace_)),
      du_(mfem::ParGridFunction(HCurlFESpace_)) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);

  true_offsets.SetSize(3);
  true_offsets[0] = 0;
  true_offsets[1] = H1FESpace_->GetVSize();
  true_offsets[2] = HCurlFESpace_->GetVSize();
  true_offsets.PartialSum();

  this->height = true_offsets[2];
  this->width = true_offsets[2];
}

void HCurlSolver::Init(mfem::Vector &X) {
  // Define material property coefficients
  dtCoef = mfem::ConstantCoefficient(1.0);
  oneCoef = mfem::ConstantCoefficient(1.0);
  SetMaterialCoefficients(_domain_properties);
  dtAlphaCoef = new mfem::TransformedCoefficient(&dtCoef, alphaCoef, prodFunc);
  SetVariableNames();

  SetSourceCoefficient(_domain_properties);
  if (sourceVecCoef) {
    buildSource();
  }

  // a0(p, p') = (β ∇ p, ∇ p')
  a0 = new mfem::ParBilinearForm(H1FESpace_);
  a0->AddDomainIntegrator(new mfem::DiffusionIntegrator(*betaCoef));
  a0->Assemble();

  this->buildM1(betaCoef);    // (βu, u')
  this->buildCurl(alphaCoef); // (α∇×u_{n}, ∇×u')
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

  p_.ProjectCoefficient(Zero);
  u_.ProjectCoefficient(Zero_vec);
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
void HCurlSolver::ImplicitSolve(const double dt, const mfem::Vector &X,
                                mfem::Vector &dX_dt) {
  dX_dt = 0.0;
  dtCoef.constant = dt;

  p_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  u_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[1]);

  dp_.MakeRef(H1FESpace_, dX_dt, true_offsets[0]);
  du_.MakeRef(HCurlFESpace_, dX_dt, true_offsets[1]);

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
  // (α∇×u_{n}, ∇×u') + (αdt∇×du/dt_{n+1}, ∇×u') + (βdu/dt_{n+1}, u')
  // - (s0_{n+1}, u') - <(α∇×u_{n+1}) × n, u'> = 0

  // a1(du/dt_{n+1}, u') = b1(u')
  // a1(u, u') = (βu, u') + (αdt∇×u, ∇×u')
  // b1(u') = (s0_{n+1}, u') - (α∇×u_{n}, ∇×u') + <(α∇×u_{n+1}) × n, u'>

  // (α∇×u_{n}, ∇×u')
  // v_ is a grid function but curlCurl is not parallel assembled so is OK
  curlCurl->MultTranspose(u_, *b1);
  *b1 *= -1.0;

  // use du_ as a temporary
  // (s0_{n+1}, u')
  grad->Mult(p_, du_);
  m1->AddMult(du_, *b1, 1.0);
  if (src_gf) {
    src_gf->ProjectCoefficient(*sourceVecCoef);
    // Compute the discretely divergence-free portion of src_gf
    divFreeProj->Mult(*src_gf, *div_free_src_gf);
    // Compute the dual of div_free_src_gf
    hCurlMass->AddMult(*div_free_src_gf, *b1);
  }

  mfem::ParGridFunction J_gf(HCurlFESpace_);
  mfem::Array<int> ess_tdof_list;
  J_gf = 0.0;
  _bc_map.applyEssentialBCs(u_name, ess_tdof_list, J_gf, pmesh_);
  _bc_map.applyIntegratedBCs(u_name, *b1, pmesh_);

  // a1(du/dt_{n+1}, u') = (βdu/dt_{n+1}, u') + (αdt∇×du/dt_{n+1}, ∇×u')
  if (a1 == NULL || fabs(dt - dt_A1) > 1.0e-12 * dt) {
    this->buildA1(betaCoef, dtAlphaCoef);
  }
  a1->FormLinearSystem(ess_tdof_list, J_gf, *b1, *A1, *X1, *B1);

  // We only need to create the solver and preconditioner once
  if (ams_a1 == NULL) {
    ams_a1 = new mfem::HypreAMS(*A1, HCurlFESpace_);
    ams_a1->SetSingularProblem();
  }
  if (pcg_a1 == NULL) {
    pcg_a1 = new mfem::HyprePCG(*A1);
    pcg_a1->SetTol(1.0e-16);
    pcg_a1->SetMaxIter(1000);
    pcg_a1->SetPrintLevel(0);
    pcg_a1->SetPreconditioner(*ams_a1);
  }
  // solve the system
  pcg_a1->Mult(*B1, *X1);

  a1->RecoverFEMSolution(*X1, *b1, du_);
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
  if (curlCurl != NULL) {
    delete curlCurl;
  }

  curlCurl = new mfem::ParBilinearForm(HCurlFESpace_);
  curlCurl->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*MuInv));
  curlCurl->Assemble();

  // no ParallelAssemble since this will be applied to GridFunctions
}

void HCurlSolver::SetVariableNames() {
  p_name = "scalar_potential";
  p_display_name = "Scalar Potential";

  u_name = "h_curl_var";
  u_display_name = "H(Curl) variable";
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

void HCurlSolver::SetSourceCoefficient(
    hephaestus::DomainProperties &domain_properties) {
  sourceVecCoef = domain_properties.vector_property_map["source"];
}

void HCurlSolver::buildSource() {
  src_gf = new mfem::ParGridFunction(HCurlFESpace_);
  div_free_src_gf = new mfem::ParGridFunction(HCurlFESpace_);
  // int irOrder = H1FESpace_->GetElementTransformation(0)->OrderW() +
  //               2 * H1FESpace_->GetOrder();
  int irOrder = H1FESpace_->GetElementTransformation(0)->OrderW() + 2 * 2;
  divFreeProj = new mfem::common::DivergenceFreeProjector(
      *H1FESpace_, *HCurlFESpace_, irOrder, NULL, NULL, NULL);
  hCurlMass = new mfem::ParBilinearForm(HCurlFESpace_);
  hCurlMass->AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
  hCurlMass->Assemble();
}

void HCurlSolver::RegisterOutputFields(mfem::DataCollection *dc_) {
  dc_->RegisterField(u_name, &u_);
  dc_->RegisterField(p_name, &p_);
  dc_->RegisterField("source", div_free_src_gf);
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

  socks_[p_name] = new mfem::socketstream;
  socks_[p_name]->precision(8);

  socks_["source"] = new mfem::socketstream;
  socks_["source"]->precision(8);

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

  mfem::common::VisualizeField(*socks_[p_name], vishost, visport, p_,
                               p_display_name.c_str(), Wx, Wy, Ww, Wh);
  Wx += offx;
  mfem::common::VisualizeField(*socks_["source"], vishost, visport,
                               *div_free_src_gf, "source", Wx, Wy, Ww, Wh);
  Wx += offx;
}

} // namespace hephaestus
