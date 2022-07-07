// Strong form
// ∇⋅v0 = 0
// ∇×(α∇×u) - βdu/dt = v0

// Weak form
// -(v0, ∇ p') + <n.v0, p'> = 0
// (α∇×u, ∇×u') - (βdu/dt, u') - (v0, u') - <(α∇×u) × n, u'> = 0

// v, v0 ∈ H(div) source field
// u ∈ H(curl)
// p ∈ H1

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
      v_(mfem::ParGridFunction(H1FESpace_)),
      e_(mfem::ParGridFunction(HCurlFESpace_)),
      b_(mfem::ParGridFunction(HDivFESpace_)),
      dv_(mfem::ParGridFunction(H1FESpace_)),
      de_(mfem::ParGridFunction(HCurlFESpace_)),
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

  // Define material property coefficients
  dtCoef = mfem::ConstantCoefficient(1.0);
  oneCoef = mfem::ConstantCoefficient(1.0);
  SetMaterialCoefficients(domain_properties);
  SetVariableNames();

  // Variables
  // v_, "electric_potential"
  // e_, "electric_field"
  // b_, "magnetic_flux_density", "Magnetic Flux Density (B)" (auxvar)

  // Coefficients
  // σ, "electrical_conductivity"
  // mu, "magnetic_permeability"

  // Bilinear for divergence free source field solve
  // -(J0, ∇ V') + <n.J, V'> = 0  where J0 = -σ∇V
  a0 = new mfem::ParBilinearForm(H1FESpace_);
  a0->AddDomainIntegrator(new mfem::DiffusionIntegrator(*sigmaCoef));
  a0->Assemble();

  // (ν∇×E, ∇×E') - (σE, E') - (J0, E') - <(ν∇×E) × n, E'> = 0
  // a1 = -σ(dE, E') + dt (ν∇×E, ∇×E')
  // b1 = dt [(J0, E') + <(ν∇×E) × n, E'>]

  this->buildM1(sigmaCoef);
  this->buildCurl(muInvCoef);
  this->buildGrad();
  b0 = new mfem::ParLinearForm(H1FESpace_);
  b1 = new mfem::ParLinearForm(HCurlFESpace_);
  A0 = new mfem::HypreParMatrix;
  A1 = new mfem::HypreParMatrix;
  X0 = new mfem::Vector;
  X1 = new mfem::Vector;
  B0 = new mfem::Vector;
  B1 = new mfem::Vector;

} // namespace hephaestus

/*
This is the main computational code that computes dX/dt implicitly
where X is the state vector containing A and V

(M1+dt S1) E = WeakCurl^T B + Grad V
        S0 V = 0
*/
void HCurlSolver::ImplicitSolve(const double dt, const mfem::Vector &X,
                                mfem::Vector &dX_dt) {
  dX_dt = 0.0;
  dtCoef.constant = dt;

  v_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  e_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[1]);
  b_.MakeRef(HDivFESpace_, const_cast<mfem::Vector &>(X), true_offsets[2]);

  dv_.MakeRef(H1FESpace_, dX_dt, true_offsets[0]);
  de_.MakeRef(HCurlFESpace_, dX_dt, true_offsets[1]);
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
  // v1 = <1/mu v, curl u> B
  // B is a grid function but weakCurl is not parallel assembled so is OK
  weakCurl->MultTranspose(b_, *b1);

  // use e_ as a temporary, E = Grad P
  // b1 = -dt * Grad V
  grad->Mult(v_, e_);
  m1->AddMult(e_, *b1, 1.0);

  mfem::ParGridFunction J_gf(HCurlFESpace_);
  mfem::Array<int> ess_tdof_list;
  J_gf = 0.0;
  _bc_map.applyEssentialBCs(u_name, ess_tdof_list, J_gf, pmesh_);
  _bc_map.applyIntegratedBCs(u_name, *b1, pmesh_);
  if (a1 == NULL || fabs(dt - dt_A1) > 1.0e-12 * dt) {
    this->buildA1(sigmaCoef, dtMuInvCoef);
  }
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

  a1->RecoverFEMSolution(*X1, *b1, e_);
  de_ = 0.0;

  // the total field is E_tot = E_ind - Grad Phi
  // so we need to subtract out Grad Phi
  // E = E - grad (P)
  // note grad maps GF to GF
  grad->AddMult(v_, e_, -1.0);

  // Compute dB/dt = -Curl(E_{n+1})
  // note curl maps GF to GF
  curl->Mult(e_, db_);
  db_ *= -1.0;
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
void HCurlSolver::Init(mfem::Vector &X) {
  mfem::Vector zero_vec(3);
  zero_vec = 0.0;
  mfem::VectorConstantCoefficient Zero_vec(zero_vec);
  mfem::ConstantCoefficient Zero(0.0);

  v_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  e_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[1]);
  b_.MakeRef(HDivFESpace_, const_cast<mfem::Vector &>(X), true_offsets[2]);

  v_.ProjectCoefficient(Zero);
  e_.ProjectCoefficient(Zero_vec);
  b_.ProjectCoefficient(Zero_vec);
}

void HCurlSolver::SetVariableNames() {
  p_name = "electric_potential";
  p_display_name = "Scalar Potential (V)";

  u_name = "electric_field";
  u_display_name = "Electric Field (E)";

  v_name = "magnetic_flux_density";
  v_display_name = "Magnetic Flux Density (B)";
}

void HCurlSolver::SetMaterialCoefficients(
    hephaestus::DomainProperties &domain_properties) {
  sigmaCoef = domain_properties.scalar_property_map["electrical_conductivity"];
  muCoef = domain_properties.scalar_property_map["magnetic_permeability"];
  muInvCoef = new mfem::TransformedCoefficient(&oneCoef, muCoef, fracFunc);
  dtMuInvCoef = new mfem::TransformedCoefficient(&dtCoef, muInvCoef, prodFunc);
}

void HCurlSolver::RegisterOutputFields(mfem::DataCollection *dc_) {
  dc_->RegisterField(u_name, &e_);
  dc_->RegisterField(v_name, &b_);
  dc_->RegisterField(p_name, &v_);
}

void HCurlSolver::WriteConsoleSummary(double t, int it) {
  // Write a summary of the timestep to console.

  // Output Ohmic losses to console
  double el = this->ElectricLosses();
  if (myid_ == 0) {
    std::cout << std::fixed;
    std::cout << "step " << std::setw(6) << it << ",\tt = " << std::setw(6)
              << std::setprecision(3) << t
              << ",\tdot(E, J) = " << std::setprecision(8) << el << std::endl;
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

  mfem::common::VisualizeField(*socks_[u_name], vishost, visport, e_,
                               u_display_name.c_str(), Wx, Wy, Ww, Wh);
  Wx += offx;

  mfem::common::VisualizeField(*socks_[v_name], vishost, visport, b_,
                               v_display_name.c_str(), Wx, Wy, Ww, Wh);
  Wx += offx;

  mfem::common::VisualizeField(*socks_[p_name], vishost, visport, v_,
                               p_display_name.c_str(), Wx, Wy, Ww, Wh);
  Wx += offx;
}

double HCurlSolver::ElectricLosses() const {
  double el = m1->InnerProduct(e_, e_);

  double global_el;
  MPI_Allreduce(&el, &global_el, 1, MPI_DOUBLE, MPI_SUM,
                m1->ParFESpace()->GetComm());

  return el;
}
} // namespace hephaestus