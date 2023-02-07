// Linear Elastic Solver
// Mostly ripped from example two
// *formulation*

#include "linear_elastic_solver.hpp"

namespace hephaestus {

LinearElasticSolver::LinearElasticSolver(
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
      a1(NULL), a1_solver(NULL),
      u_(mfem::ParGridFunction(H1FESpace_)),
      du_(mfem::ParGridFunction(H1FESpace_)) {
  // Initialize MPI variables
  MPI_Comm_size(pmesh.GetComm(), &num_procs_);
  MPI_Comm_rank(pmesh.GetComm(), &myid_);

  true_offsets.SetSize(2);
  true_offsets[0] = 0;
  true_offsets[1] = H1FESpace_->GetVSize();
  true_offsets.PartialSum();

  this->height = true_offsets[1];
  this->width = true_offsets[1];

  HYPRE_BigInt size_h1 = H1FESpace_->GlobalTrueVSize();
  if (myid_ == 0) {
    std::cout << "Total number of         DOFs: " << size_h1 << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "Total number of H1      DOFs: " << size_h1 << std::endl;
    std::cout << "------------------------------------" << std::endl;
  }
}

void LinearElasticSolver::Init(mfem::Vector &X) {  
  //Register the pargridfunction associated with strain in the variables NamedFieldsMap
  RegisterVariables();
  _fespaces.Register("_H1FESpace", H1FESpace_, false);

  // Initialise coefficients and store in DomainProperties
  SetMaterialCoefficients(_domain_properties);

  //Allocate meory of bilinear forms, linear forms and matrices

  // Set up pull force for linear elastic problem
  mfem::VectorArrayCoefficient f(pmesh_->Dimension());
  // Set all components to 0
  for (int i = 0; i < pmesh_->Dimension()-1; i++)
  {
    f.Set(i, new mfem::ConstantCoefficient(0.0));
  }
  
  // Set Z component
  mfem::Vector pull_force(pmesh_->bdr_attributes.Max());
  pull_force = 0.0;
  pull_force(1) = -1.0e-2;
  f.Set(pmesh_->Dimension()-1, new mfem::PWConstCoefficient(pull_force));
  
  b1 = new mfem::ParLinearForm(H1FESpace_);
  A1 = new mfem::HypreParMatrix;
  X1 = new mfem::Vector;
  B1 = new mfem::Vector;

  this->buildA1();
  // Set up linear form for linear elastic problem
  b1->AddBoundaryIntegrator(new mfem::VectorBoundaryLFIntegrator(f));
  b1->Assemble();

  u_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);

  //Initial state
  mfem::Vector zero_vec;   
  zero_vec = 0.0;
  mfem::VectorConstantCoefficient Zero_vec(zero_vec);
  u_.ProjectCoefficient(Zero_vec);
}

/*
This is the main computational code that computes dX/dt implicitly
where X is the state vector containing p, u and v.
*/


// Mult method is called instead of implicitSolve for steady state solves
void LinearElasticSolver::Mult(const mfem::Vector &X)
{
  //Update gridfunctions to ref solution vector X and time derivative time dX/dt 
  u_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);

  // Update linear form - can be given from bcmap 
  *b1 = 0.0;
  _sources.ApplyKernels(b1);
  mfem::ParGridFunction x_gf(H1FESpace_);
  mfem::Array<int> ess_tdof_list, ess_bdr(pmesh_->bdr_attributes.Max());
  ess_bdr = 0;
  ess_bdr[0] = 1;
  H1FESpace_->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  x_gf = 0.0;

   //Form linear system from blf and lf
  a1->FormLinearSystem(ess_tdof_list, x_gf, *b1, *A1, *X1, *B1);

  // We only need to create the solver and preconditioner once
  mfem::HypreBoomerAMG *amg = new mfem::HypreBoomerAMG(*A1);
  
  amg->SetSystemsOptions(pmesh_->Dimension(), false);
  mfem::HyprePCG *pcg = new mfem::HyprePCG(*A1);
  pcg->SetTol(1e-8);
  pcg->SetMaxIter(500);
  pcg->SetPrintLevel(2);
  pcg->SetPreconditioner(*amg);
  pcg->Mult(*B1, *X1);

  a1->RecoverFEMSolution(*X1, *b1, du_);

  pmesh_->SetNodalFESpace(H1FESpace_);
}


void LinearElasticSolver::ImplicitSolve(const double dt, const mfem::Vector &X,
                                mfem::Vector &dX_dt) {
  dX_dt = 0.0;
  dtCoef.constant = dt;

  //Update gridfunctions to ref solution vector X and time derivative time dX/dt 
  u_.MakeRef(H1FESpace_, const_cast<mfem::Vector &>(X), true_offsets[0]);
  // du_.MakeRef(H1FESpace_, dX_dt, true_offsets[0]);

  // Commented out as this is not currently time dependent 
  // _domain_properties.SetTime(this->GetTime());

  // Update linear form - can be given from bcmap 
  *b1 = 0.0;
  _sources.ApplyKernels(b1);
  mfem::ParGridFunction x_gf(H1FESpace_);
  mfem::Array<int> ess_tdof_list, ess_bdr(pmesh_->bdr_attributes.Max());
  ess_bdr = 0;
  ess_bdr[0] = 1;
  H1FESpace_->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  x_gf = 0.0;

   //Form linear system from blf and lf
  a1->FormLinearSystem(ess_tdof_list, x_gf, *b1, *A1, *X1, *B1);

  // We only need to create the solver and preconditioner once

  mfem::HypreBoomerAMG *amg = new mfem::HypreBoomerAMG(*A1);
  // if (amg_elast && !a1->StaticCondensationIsEnabled())
  // {
  //   amg->SetElasticityOptions(H1FESpace_);
  // }
  // else
  // {
    amg->SetSystemsOptions(pmesh_->Dimension(), false);
  // }
  mfem::HyprePCG *pcg = new mfem::HyprePCG(*A1);
  pcg->SetTol(1e-8);
  pcg->SetMaxIter(500);
  pcg->SetPrintLevel(2);
  pcg->SetPreconditioner(*amg);
  pcg->Mult(*B1, *X1);

  a1->RecoverFEMSolution(*X1, *b1, du_);

  //Deform mesh?
  //  if (!use_nodal_fespace)
  //  {
      pmesh_->SetNodalFESpace(H1FESpace_);
  //  }

}

void LinearElasticSolver::buildA1() {
  if (a1 != NULL) {
    delete a1;
  }

  // First create and assemble the bilinear form.  For now we assume the mesh
  // isn't moving, the materials are time independent, and dt is constant. So
  // we only need to do this once.
  a1 = new mfem::ParBilinearForm(H1FESpace_);
  a1->AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda_func, mu_func));
  a1->Assemble();

  // Don't finalize or parallel assemble this is done in FormLinearSystem.
  dt_A1 = dtCoef.constant;
}

// Register strain variable "u_" in _variables
void LinearElasticSolver::RegisterVariables() {
  u_name = "strain";
  _variables.Register(u_name, &u_, false);
}

void LinearElasticSolver::SetMaterialCoefficients(
    hephaestus::DomainProperties &domain_properties) {

  // Set up lame coefficients
  mfem::Vector lambda(pmesh_->attributes.Max());
  std::cout << "Number of attrs: " << pmesh_->attributes.Max() << std::endl;
  lambda = 1.0;
  lambda(0) = lambda(1)*50;
  lambda_func = mfem::PWConstCoefficient(lambda);

  mfem::Vector mu(pmesh_->attributes.Max());
  mu = 1.0;
  mu(0) = mu(1)*50;
  mu_func = mfem::PWConstCoefficient(mu);
}

// Put our strain data into datacollection object dc_
void LinearElasticSolver::RegisterOutputFields(mfem::DataCollection *dc_) {
  dc_->SetMesh(pmesh_);
  for (auto var = _variables.begin(); var != _variables.end(); ++var) {
    dc_->RegisterField(var->first, var->second);
  }
}

void LinearElasticSolver::WriteConsoleSummary(double t, int it) {
  // Write a summary of the timestep to console.
  if (myid_ == 0) {
    std::cout << std::fixed;
    std::cout << "step " << std::setw(6) << it << ",\tt = " << std::setw(6)
              << std::setprecision(3) << t << std::endl;
  }
}

void LinearElasticSolver::WriteOutputFields(mfem::DataCollection *dc_, int it) {
  if (dc_) {
    dc_->SetCycle(it);
    dc_->SetTime(t);
    dc_->Save();
  }
}

void LinearElasticSolver::InitializeGLVis() {
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

void LinearElasticSolver::DisplayToGLVis() {
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
