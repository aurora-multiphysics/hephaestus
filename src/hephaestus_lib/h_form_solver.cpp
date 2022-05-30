#include "h_form_solver.hpp"

namespace hephaestus {
void hdot_bc_func(const mfem::Vector &x, mfem::Vector &E) { E = 0.0; }

HFormSolver::HFormSolver(mfem::ParMesh &pmesh, int order,
                         hephaestus::BCMap &bc_map,
                         hephaestus::DomainProperties &domain_properties)
    : myid_(0), num_procs_(1), pmesh_(&pmesh), _bc_map(bc_map),
      _domain_properties(domain_properties),
      HCurlFESpace_(
          new mfem::common::ND_ParFESpace(&pmesh, order, pmesh.Dimension())),
      amg_a0(NULL), pcg_a0(NULL), ams_a1(NULL), pcg_a1(NULL),
      h_(mfem::ParGridFunction(HCurlFESpace_)),
      dh_(mfem::ParGridFunction(HCurlFESpace_)) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);
  this->height = HCurlFESpace_->GlobalTrueVSize();
  this->width = HCurlFESpace_->GlobalTrueVSize();

  // Define material property coefficients
  dt = 0.5;
  rho = 1.0;
  muCoef = domain_properties.getGlobalScalarProperty(
      std::string("magnetic_permeability"));
  dtRhoCoef = mfem::ConstantCoefficient(dt * rho);

  // (ρ∇×H, ∇×H') + μ(dH, H') - (dB0/dt, H') - <(ρ∇×H) × n, H'> = 0
  // a1 = μ(dH, H') + dt (ρ∇×H, ∇×H')
  // b1 = (B0, H') + dt<(ρ∇×H) × n, H'>
  a0 = new mfem::ParBilinearForm(HCurlFESpace_);
  a0->AddDomainIntegrator(new mfem::CurlCurlIntegrator(dtRhoCoef));
  a0->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(muCoef));
  a0->Assemble();

  this->buildM1(muCoef);
  this->buildGrad();
  b0 = new mfem::ParGridFunction(HCurlFESpace_);
  A0 = new mfem::HypreParMatrix;
  X0 = new mfem::Vector;
  B0 = new mfem::Vector;
}

/*
This is the main computational code that computes dX/dt implicitly
where X is the state vector containing A and V

(M1+dt S1) A = WeakCurl^T B + Grad V
        S0 V = 0
*/
void HFormSolver::ImplicitSolve(const double dt, const mfem::Vector &X,
                                mfem::Vector &dX_dt) {
  dX_dt = 0.0;

  int Vsize_h1 = H1FESpace_->GetVSize();
  int Vsize_nd = HCurlFESpace_->GetVSize();

  mfem::Array<int> true_offset(3);
  true_offset[0] = 0;
  true_offset[1] = true_offset[0] + Vsize_nd;

  mfem::Vector *xptr = (mfem::Vector *)&X;
  h_.MakeRef(HCurlFESpace_, *xptr, true_offset[0]);
  dh_.MakeRef(HCurlFESpace_, dX_dt, true_offset[0]);

  // Dirichlet conditions
  mfem::Array<int> ess_bdr = _bc_map["tangential_dEdt"]->getMarkers(*pmesh_);
  mfem::VectorFunctionCoefficient Hdot_bc(3, hdot_bc);
  mfem::ParGridFunction J_gf(HCurlFESpace_);
  J_gf = 0.0;
  J_gf.ProjectBdrCoefficientTangent(Hdot_bc, ess_bdr);

  // form the linear system, including eliminating essential BC's and
  // applying
  // static condensation. The system to solve is A1 X1 = B1
  mfem::Array<int> ess_tdof_list;
  HCurlFESpace_->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

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

  a1->RecoverFEMSolution(*X1, *b1, h_);
}

void HFormSolver::buildM1(mfem::PWCoefficient &Sigma) {
  if (m1 != NULL) {
    delete m1;
  }

  m1 = new mfem::ParBilinearForm(HCurlFESpace_);
  m1->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(Sigma));
  m1->Assemble();

  // Don't finalize or parallel assemble this is done in FormLinearSystem.
}

void HFormSolver::buildGrad() {
  if (grad != NULL) {
    delete grad;
  }

  grad = new mfem::ParDiscreteLinearOperator(H1FESpace_, HCurlFESpace_);
  grad->AddDomainInterpolator(new mfem::GradientInterpolator());
  grad->Assemble();

  // no ParallelAssemble since this will be applied to GridFunctions
}

void HFormSolver::Init(mfem::Vector &X) {
  mfem::Vector zero_vec(3);
  zero_vec = 0.0;
  mfem::VectorConstantCoefficient Zero_vec(zero_vec);
  mfem::ConstantCoefficient Zero(0.0);

  int Vsize_h1 = H1FESpace_->GetVSize();
  int Vsize_nd = HCurlFESpace_->GetVSize();

  mfem::Array<int> true_offset(3);
  true_offset[0] = 0;
  true_offset[1] = true_offset[0] + Vsize_h1;
  true_offset[2] = true_offset[1] + Vsize_nd;

  mfem::Vector *xptr = (mfem::Vector *)&X;
  v_.MakeRef(H1FESpace_, *xptr, true_offset[0]);
  h_.MakeRef(HCurlFESpace_, *xptr, true_offset[1]);

  v_.ProjectCoefficient(Zero);
  h_.ProjectCoefficient(Zero_vec);
}

void HFormSolver::RegisterVisItFields(mfem::VisItDataCollection &visit_dc) {
  visit_dc_ = &visit_dc;

  visit_dc.RegisterField("A", &h_);
  visit_dc.RegisterField("V", &v_);
}

void HFormSolver::WriteVisItFields(int it) {
  if (visit_dc_) {
    if (myid_ == 0) {
      std::cout << "Writing VisIt files ..." << std::flush;
    }

    visit_dc_->SetCycle(it);
    visit_dc_->SetTime(t);
    visit_dc_->Save();

    if (myid_ == 0) {
      std::cout << " done." << std::endl;
    }
  }
}

void HFormSolver::InitializeGLVis() {
  if (myid_ == 0) {
    std::cout << "Opening GLVis sockets." << std::endl;
  }

  socks_["A"] = new mfem::socketstream;
  socks_["A"]->precision(8);

  socks_["V"] = new mfem::socketstream;
  socks_["V"]->precision(8);

  if (myid_ == 0) {
    std::cout << "GLVis sockets open." << std::endl;
  }
}

void HFormSolver::DisplayToGLVis() {
  if (myid_ == 0) {
    std::cout << "Sending data to GLVis ..." << std::flush;
  }

  char vishost[] = "localhost";
  int visport = 19916;

  int Wx = 0, Wy = 0;                 // window position
  int Ww = 350, Wh = 350;             // window size
  int offx = Ww + 10, offy = Wh + 45; // window offsets

  mfem::common::VisualizeField(*socks_["A"], vishost, visport, h_,
                               "Vector Potential (A)", Wx, Wy, Ww, Wh);
  Wx += offx;

  mfem::common::VisualizeField(*socks_["V"], vishost, visport, v_,
                               "Scalar Potential (V)", Wx, Wy, Ww, Wh);
  Wx += offx;

  if (myid_ == 0) {
    std::cout << " done." << std::endl;
  }
}

} // namespace hephaestus