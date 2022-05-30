#include "e_solver.hpp"

namespace hephaestus {
void edot_bc(const mfem::Vector &x, mfem::Vector &E) { E = 0.0; }

ESolver::ESolver(mfem::ParMesh &pmesh, int order, hephaestus::BCMap &bc_map,
                 hephaestus::DomainProperties &domain_properties)
    : myid_(0), num_procs_(1), pmesh_(&pmesh), _bc_map(bc_map),
      _domain_properties(domain_properties),
      H1FESpace_(
          new mfem::common::H1_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HCurlFESpace_(
          new mfem::common::ND_ParFESpace(&pmesh, order, pmesh.Dimension())),
      amg_a0(NULL), pcg_a0(NULL), ams_a1(NULL), pcg_a1(NULL), m1(NULL),
      grad(NULL), v_(mfem::ParGridFunction(H1FESpace_)),
      e_(mfem::ParGridFunction(HCurlFESpace_)),
      dv_(mfem::ParGridFunction(H1FESpace_)),
      de_(mfem::ParGridFunction(HCurlFESpace_)) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);
  this->height =
      HCurlFESpace_->GlobalTrueVSize() + H1FESpace_->GlobalTrueVSize();
  this->width =
      HCurlFESpace_->GlobalTrueVSize() + H1FESpace_->GlobalTrueVSize();
  // Define compatible parallel finite element spaces on the parallel
  // mesh. Here we use arbitrary order H1, Nedelec, and Raviart-Thomas finite
  // elements.
  // H1FESpace_ =
  //     new mfem::common::H1_ParFESpace(&pmesh, order, pmesh.Dimension());
  // HCurlFESpace_ =
  //     new mfem::common::ND_ParFESpace(&pmesh, order, pmesh.Dimension());
  // HDivFESpace_ =
  //     new mfem::common::RT_ParFESpace(&pmesh, order, pmesh.Dimension());

  // Define gridfunctions
  // v_ = new mfem::ParGridFunction(H1FESpace_);     // Scalar Potential (H1)
  // e_ = new mfem::ParGridFunction(HCurlFESpace_);  // Vector Potential (HCurl)
  // dv_ = new mfem::ParGridFunction(H1FESpace_);    // Scalar Potential (H1)
  // de_ = new mfem::ParGridFunction(HCurlFESpace_); // Vector Potential (HCurl)

  // e_ = new mfem::ParGridFunction(HCurlFESpace_); // Electric Field (HCurl)
  // b_ = new mfem::ParGridFunction(HDivFESpace_);  // Magnetic Flux (HDiv)

  // Define material property coefficients
  dt = 0.5;
  mu = 1.0;
  sigmaCoef = domain_properties.getGlobalScalarProperty(
      std::string("electrical_conductivity"));
  dtMuInvCoef = mfem::ConstantCoefficient(dt * mu);

  // (σ(dA/dt + ∇ V), ∇ V') + <n.J, V'> = 0
  a0 = new mfem::ParBilinearForm(H1FESpace_);
  a0->AddDomainIntegrator(new mfem::DiffusionIntegrator(sigmaCoef));
  a0->Assemble();

  // σ[(e_{n+1}-e_n)/dt] = f(t_{n+1}, e_{n+1})
  // σ(dA/dt) = σ((e_{n+1}-e_n)/dt)  --> σ e_{n+1} = σ dt e_n
  // (ν∇×A, ∇×A') + (σ(dA/dt + ∇ V), A') - (J0, A') - <(ν∇×A) × n, A'> = 0
  // a1 = σ(dA, A') + dt (ν∇×A, ∇×A')
  // b1 = dt [-(σ ∇ V, A') + (J0, A') + <(ν∇×A) × n, A'>]
  a1 = new mfem::ParBilinearForm(HCurlFESpace_);
  a1->AddDomainIntegrator(new mfem::CurlCurlIntegrator(dtMuInvCoef));
  a1->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(sigmaCoef));
  a1->Assemble();

  this->buildM1(sigmaCoef);
  this->buildGrad();
  b0 = new mfem::ParGridFunction(H1FESpace_);
  b1 = new mfem::ParGridFunction(HCurlFESpace_);
  A0 = new mfem::HypreParMatrix;
  A1 = new mfem::HypreParMatrix;
  X0 = new mfem::Vector;
  X1 = new mfem::Vector;
  B0 = new mfem::Vector;
  B1 = new mfem::Vector;

  // // (ν∇×A, ∇×A') + (σ(dA/dt + ∇ V), A') - (J0, A') - <(ν∇×A) × n, A'> = 0
  // // (ν∇×A, ∇×A')

  // // (ν∇×A, ∇×A') + (σ(dA/dt + ∇ V), A') - (J0, A') - <(ν∇×A) × n, A'> = 0
  // // (ν∇×A, ∇×A')
  // curlMuInvCurl_ = new mfem::ParBilinearForm(HCurlFESpace_);
  // curlMuInvCurl_->AddDomainIntegrator(new
  // mfem::CurlCurlIntegrator(muInvCoef));

  // // (σ dA/dt, A')
  // hCurlMass_ = new mfem::ParBilinearForm(HCurlFESpace_);
  // hCurlMass_->AddDomainIntegrator(new
  // mfem::VectorFEMassIntegrator(sigmaCoef));

  // // (σ ∇ V, A')
  // sigmaGradH1HCurl_ = new mfem::ParMixedBilinearForm(H1FESpace_,
  // HCurlFESpace_); sigmaGradH1HCurl_->AddDomainIntegrator(
  //     new mfem::MixedVectorGradientIntegrator(sigmaCoef));

  // (σ(dA/dt + ∇ V), ∇ V') + <n.J, V'> = 0
  // BilinearFormIntegrator *hCurlMassInteg = new VectorFEMassIntegrator;
  // hCurlMass_ = new mfem::ParBilinearForm(HCurlFESpace_);
  // hCurlMass_->AddDomainIntegrator(hCurlMassInteg);

  // (σ(dA/dt + ∇ V), ∇ V') + <n.J, V'> = 0

  // this->buildA0(*sigma);
  // this->buildM1(*sigma);
  // this->buildS1(1.0 / mu);
  // this->buildCurl(1.0 / mu);
  // this->buildGrad();
} // namespace hephaestus

/*
This is the main computational code that computes dX/dt implicitly
where X is the state vector containing A and V

(M1+dt S1) A = WeakCurl^T B + Grad V
        S0 V = 0
*/
void ESolver::ImplicitSolve(const double dt, const mfem::Vector &X,
                            mfem::Vector &dX_dt) {
  dX_dt = 0.0;

  int Vsize_h1 = H1FESpace_->GetVSize();
  int Vsize_nd = HCurlFESpace_->GetVSize();

  mfem::Array<int> true_offset(3);
  true_offset[0] = 0;
  true_offset[1] = true_offset[0] + Vsize_h1;
  true_offset[2] = true_offset[1] + Vsize_nd;

  mfem::Vector *xptr = (mfem::Vector *)&X;
  v_.MakeRef(H1FESpace_, *xptr, true_offset[0]);
  e_.MakeRef(HCurlFESpace_, *xptr, true_offset[1]);

  dv_.MakeRef(H1FESpace_, dX_dt, true_offset[0]);
  de_.MakeRef(HCurlFESpace_, dX_dt, true_offset[1]);

  // form the Laplacian and solve it
  mfem::ParGridFunction Phi_gf(H1FESpace_);

  // p_bc is given function defining electrostatic potential on surface
  mfem::Array<int> poisson_ess_bdr =
      _bc_map["electric_potential"]->getMarkers(*pmesh_);

  hephaestus::FunctionDirichletBC *potential_bc =
      dynamic_cast<hephaestus::FunctionDirichletBC *>(
          _bc_map["electric_potential"]);
  mfem::FunctionCoefficient voltage = *potential_bc->coeff;

  voltage.SetTime(this->GetTime());
  Phi_gf = 0.0;

  // the function below is currently not fully supported on AMR meshes
  // Phi_gf.ProjectBdrCoefficient(voltage,poisson_ess_bdr);

  // this is a hack to get around the above issue
  Phi_gf.ProjectCoefficient(voltage);
  // end of hack

  // apply essential BC's and apply static condensation, the new system to
  // solve is A0 X0 = B0
  mfem::Array<int> poisson_ess_tdof_list;
  H1FESpace_->GetEssentialTrueDofs(poisson_ess_bdr, poisson_ess_tdof_list);

  *b0 = 0.0;
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
  // use e_ as a temporary, E = Grad P
  // b1 = -dt * sigma * Grad V
  grad->Mult(v_, e_);
  m1->AddMult(e_, *b1, -dt);

  mfem::Array<int> ess_bdr = _bc_map["tangential_dEdt"]->getMarkers(*pmesh_);
  mfem::VectorFunctionCoefficient Jdot(3, edot_bc);
  mfem::ParGridFunction J_gf(HCurlFESpace_);
  J_gf = 0.0;
  J_gf.ProjectBdrCoefficientTangent(Jdot, ess_bdr);

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

  a1->RecoverFEMSolution(*X1, *b1, e_);
}

void ESolver::buildM1(mfem::PWCoefficient &Sigma) {
  if (m1 != NULL) {
    delete m1;
  }

  m1 = new mfem::ParBilinearForm(HCurlFESpace_);
  m1->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(Sigma));
  m1->Assemble();

  // Don't finalize or parallel assemble this is done in FormLinearSystem.
}

void ESolver::buildGrad() {
  if (grad != NULL) {
    delete grad;
  }

  grad = new mfem::ParDiscreteLinearOperator(H1FESpace_, HCurlFESpace_);
  grad->AddDomainInterpolator(new mfem::GradientInterpolator());
  grad->Assemble();

  // no ParallelAssemble since this will be applied to GridFunctions
}

void ESolver::Init(mfem::Vector &X) {
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
  e_.MakeRef(HCurlFESpace_, *xptr, true_offset[1]);

  v_.ProjectCoefficient(Zero);
  e_.ProjectCoefficient(Zero_vec);
}

void ESolver::RegisterVisItFields(mfem::VisItDataCollection &visit_dc) {
  visit_dc_ = &visit_dc;

  visit_dc.RegisterField("A", &e_);
  visit_dc.RegisterField("V", &v_);
}

void ESolver::WriteVisItFields(int it) {
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

void ESolver::InitializeGLVis() {
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

void ESolver::DisplayToGLVis() {
  if (myid_ == 0) {
    std::cout << "Sending data to GLVis ..." << std::flush;
  }

  char vishost[] = "localhost";
  int visport = 19916;

  int Wx = 0, Wy = 0;                 // window position
  int Ww = 350, Wh = 350;             // window size
  int offx = Ww + 10, offy = Wh + 45; // window offsets

  mfem::common::VisualizeField(*socks_["A"], vishost, visport, e_,
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
