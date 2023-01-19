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

HCurlSolver::HCurlSolver(
    mfem::ParMesh &pmesh, int order,
    mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> &fespaces,
    mfem::NamedFieldsMap<mfem::ParGridFunction> &variables,
    hephaestus::BCMap &bc_map, hephaestus::DomainProperties &domain_properties,
    hephaestus::Sources &sources, hephaestus::InputParameters &solver_options)
    : myid_(0), num_procs_(1), pmesh_(&pmesh), _fespaces(fespaces),
      _variables(variables), _bc_map(bc_map), _sources(sources),
      _domain_properties(domain_properties), _solver_options(solver_options),
      H1FESpace_(
          new mfem::common::H1_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HCurlFESpace_(
          new mfem::common::ND_ParFESpace(&pmesh, order, pmesh.Dimension())),
      HDivFESpace_(
          new mfem::common::RT_ParFESpace(&pmesh, order, pmesh.Dimension())),
      a1(NULL), a1_solver(NULL), curl(NULL), curlCurl(NULL),
      u_(mfem::ParGridFunction(HCurlFESpace_)),
      du_(mfem::ParGridFunction(HCurlFESpace_)),
      curl_u_(mfem::ParGridFunction(HDivFESpace_)) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);

  true_offsets.SetSize(2);
  true_offsets[0] = 0;
  true_offsets[1] = HCurlFESpace_->GetVSize();
  true_offsets.PartialSum();

  this->height = true_offsets[1];
  this->width = true_offsets[1];

  HYPRE_BigInt size_h1 = H1FESpace_->GlobalTrueVSize();
  HYPRE_BigInt size_nd = HCurlFESpace_->GlobalTrueVSize();
  if (myid_ == 0) {
    std::cout << "Total number of         DOFs: " << size_h1 + size_nd
              << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "Total number of H1      DOFs: " << size_h1 << std::endl;
    std::cout << "Total number of H(Curl) DOFs: " << size_nd << std::endl;
    std::cout << "------------------------------------" << std::endl;
  }
}

void HCurlSolver::Init(mfem::Vector &X) {
  RegisterVariables();
  _fespaces.Register("_H1FESpace", H1FESpace_, false);
  _fespaces.Register("_HCurlFESpace", HCurlFESpace_, false);
  _fespaces.Register("_HDivFESpace", HDivFESpace_, false);

  // Define material property coefficients
  dtCoef = mfem::ConstantCoefficient(1.0);
  oneCoef = mfem::ConstantCoefficient(1.0);
  SetMaterialCoefficients(_domain_properties);
  dtAlphaCoef = new mfem::TransformedCoefficient(&dtCoef, alphaCoef, prodFunc);

  _sources.Init(_variables, _fespaces, _bc_map, _domain_properties);

  // a0(p, p') = (β ∇ p, ∇ p')
  this->buildCurl(alphaCoef); // (α∇×u_{n}, ∇×u')
  b1 = new mfem::ParLinearForm(HCurlFESpace_);
  A1 = new mfem::HypreParMatrix;
  X1 = new mfem::Vector;
  B1 = new mfem::Vector;

  mfem::Vector zero_vec(3);
  zero_vec = 0.0;
  mfem::VectorConstantCoefficient Zero_vec(zero_vec);
  mfem::ConstantCoefficient Zero(0.0);

  u_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
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

  u_.MakeRef(HCurlFESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  du_.MakeRef(HCurlFESpace_, dX_dt, true_offsets[0]);

  _domain_properties.SetTime(this->GetTime());

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

  _sources.ApplyKernels(b1);

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
  if (a1_solver == NULL) {
    a1_solver = new hephaestus::DefaultHCurlPCGSolver(_solver_options, *A1,
                                                      HCurlFESpace_);
  }
  a1_solver->Mult(*B1, *X1);

  a1->RecoverFEMSolution(*X1, *b1, du_);

  curl->Mult(u_, curl_u_);
  curl->AddMult(du_, curl_u_, dt);
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

void HCurlSolver::buildCurl(mfem::Coefficient *MuInv) {
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

  // no ParallelAssemble since this will be applied to GridFunctions
}

void HCurlSolver::RegisterVariables() {
  u_name = "h_curl_var";
  u_display_name = "H(Curl) variable";

  _variables.Register(u_name, &u_, false);
  _variables.Register("curl h_curl_var", &curl_u_, false);
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
  dc_->SetMesh(pmesh_);
  for (auto var = _variables.begin(); var != _variables.end(); ++var) {
    dc_->RegisterField(var->first, var->second);
  }
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

  for (auto var = _variables.begin(); var != _variables.end(); ++var) {
    socks_[var->first] = new mfem::socketstream;
    socks_[var->first]->precision(8);
  }

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

  for (auto var = _variables.begin(); var != _variables.end(); ++var) {
    mfem::common::VisualizeField(*socks_[var->first], vishost, visport,
                                 *(var->second), (var->first).c_str(), Wx, Wy,
                                 Ww, Wh);
    Wx += offx;
  }
}

} // namespace hephaestus
