// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
// mpirun -np 4 valgrind ./magnetodynamic_newton -m bulk1_fine.msh

#include "hephaestus_h_form.hpp"
//#ifndef MFEM_USE_PETSC
//#error This example requires that MFEM is built with MFEM_USE_PETSC=YES
//#endif

using namespace std;
using namespace mfem;

class Hform_Operator;
class SuperconductorEJIntegrator;
class SurfaceFluxLFIntegrator;
// class PreconditionerFactory;

class ParDiscreteInterpolationOperator : public ParDiscreteLinearOperator {
public:
  ParDiscreteInterpolationOperator(ParFiniteElementSpace *dfes,
                                   ParFiniteElementSpace *rfes)
      : ParDiscreteLinearOperator(dfes, rfes) {}
  virtual ~ParDiscreteInterpolationOperator();
};
ParDiscreteInterpolationOperator::~ParDiscreteInterpolationOperator() {}

class ParDiscreteCurlOperator : public ParDiscreteInterpolationOperator {
public:
  ParDiscreteCurlOperator(ParFiniteElementSpace *dfes,
                          ParFiniteElementSpace *rfes);
};
ParDiscreteCurlOperator::ParDiscreteCurlOperator(ParFiniteElementSpace *dfes,
                                                 ParFiniteElementSpace *rfes)
    : ParDiscreteInterpolationOperator(dfes, rfes) {
  this->AddDomainInterpolator(new CurlInterpolator);
}

class MagnetodynamicSolver : public TimeDependentOperator {
protected:
  VisItDataCollection *visit_dc_;

  ParFiniteElementSpace *HCurlFESpace, *HDivFESpace, *L2FESpace;
  ParBilinearForm *Hform_LHS_linear, *Hform_massH, *General_mass;
  ParNonlinearForm *Hform_LHS_nonlinear;
  ParDiscreteCurlOperator *curl_;
  HypreParVector *H_t; // Current value of the magnetic field DoFs
  ParGridFunction *H_t1, *H_t2, *J_t1; // H field at time step 1, step 2
  Vector rho_L, rho_NL, I_direction;
  VectorConstantCoefficient *Current_direction;
  NewtonSolver newton_solver;
  Solver *J_solver;
  Solver *J_prec;
  Hform_Operator *Hform_oper;
  Coefficient *muCoef_;                 //  permeability Coefficient
  VectorFunctionCoefficient *H_BCCoef_; // Vector Potential BC Function
  PWConstCoefficient *rhoCoef_L, *rhoCoef_NL;

  int myid;
  double (*muInv)(const Vector &);
  double (*ItrFunc)(double);
  void (*H_BCFunc)(const Vector &, double, Vector &);
  double (*muInv_)(const Vector &);
  double (*ItrFunc_)(double);
  void (*H_BC)(const Vector &, double, Vector &);
  Array<int> ess_bdr; // boundary marker for outer boundaries
  Array<int> ess_domain;
  Array<int> *ess_tdof_list;
  Array<int> ess_bdr_cs; // boundary marker for conductor cross sections
public:
  MagnetodynamicSolver(ParFiniteElementSpace &HCurlFESpace_,
                       ParFiniteElementSpace &HDivFESpace_,
                       ParFiniteElementSpace &L2FESpace_, Array<int> &ess_bdr_,
                       Array<int> &ess_domain_, double (*muInv)(const Vector &),
                       double (*ItrFunc)(double),
                       void (*H_BCFunc)(const Vector &, double, Vector &));

  HYPRE_Int GetProblemSize();
  void PrintSizes();
  Vector &GetHfield() {
    return *H_t;
  } // define the field H_ to be returned from ODE_solver->step()
  void SetInitial_Hfield();
  virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);
  void SyncGridFuncs();
  void RegisterVisItFields(VisItDataCollection &visit_dc);
  void WriteVisItFields(int it = 0);
  virtual ~MagnetodynamicSolver();
};

// parameter Function

double muInv(const Vector &x) {
  const double mu0_ = 1.2566e-6;
  return mu0_;
}

// time-dependent boundary condition function on H

double ramp_current1(double t) {
  double Itr = 10 * t; // current is ramped up with time from zero 10A/s
  return Itr;
}
double ItrFunc(double t) { ramp_current1(t); }
// Time-dependent boundary condition, dH/dt, derived from known time-dependent
// background field H(t) value
void H_BCFunc(const Vector &x, double t, Vector &H_bc) {
  const double mu0_ = 1.2566e-6;
  // H_bc = 0.0; // can be changed to time dependent, or, dH/
  H_bc.SetSize(3);
  H_bc(0) = 0;
  H_bc(1) = 0;
  if (t < 0.6) {
    H_bc(2) = t * 8.5 / mu0_;
  } else if (t >= 0.6 && t < 1.1) {
    H_bc(2) = (0.5 - (t - 0.5) * 0.9) / mu0_;
  } else {
    H_bc(2) = 0.00001 / mu0_;
  }
}

int h_form_solve(int argc, char *argv[]) {
  // 1. Initialize MPI.
  int num_procs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  out << "MPI size " << num_procs << endl;
  // Parse command-line options.
  const char *mesh_file = "../data/bulk3_fine.msh";
  int Order = 1;
  int serial_ref_levels = 1;
  int parallel_ref_levels = 0;
  bool visit = true;
  double dt = 1.0e-1;
  double dtsf = 0.95;
  double ti = 0.0;
  double ts = 0.5;
  double tf = 0.1;
  const char *petscrc_file = "";

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  // args.AddOption(&petscrc_file, "-petscopts", "--petscopts",
  //   "PetscOptions file to use.");

  args.Parse();
  if (!args.Good()) {
    if (myid == 0) {
      args.PrintUsage(cout);
    }
    return 1;
  }
  if (myid == 0) {
    args.PrintOptions(cout);
  }
  // initilize PETSC
  // MFEMInitializePetsc(NULL,NULL,petscrc_file,NULL);

  // Read the (serial) mesh from the given mesh file on all processors.  We can
  // handle triangular, quadrilateral, tetrahedral, hexahedral, surface and
  // volume meshes with the same code.

  Mesh *mesh = new Mesh(mesh_file, 1, 1);

  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;

  FiniteElementCollection *fec_ND;
  FiniteElementCollection *fec_RT;
  FiniteElementCollection *fec_L2;
  fec_ND = new ND_FECollection(Order, pmesh->Dimension());
  fec_RT = new RT_FECollection(Order, pmesh->Dimension());
  fec_L2 = new L2_FECollection(Order, pmesh->Dimension());
  ParFiniteElementSpace HCurlFESpace(pmesh, fec_ND);
  ParFiniteElementSpace HDivFESpace(pmesh, fec_RT);
  ParFiniteElementSpace L2FESpace(pmesh, fec_L2);
  Vector H1(HCurlFESpace.TrueVSize());
  H1 = 0.0;
  Array<int> ess_bdr_(
      pmesh->bdr_attributes.Max()); // all boundary attributes in mesh
  Array<int> ess_domain_(
      pmesh->attributes.Max()); // all domain attributes in mesh

  int size_Hcurl = HCurlFESpace.GetTrueVSize();
  if (myid == 0) {
    cout << "Number of local H(Curl) unknowns: " << size_Hcurl << endl;
  }
  HYPRE_Int size_rt = HDivFESpace.GlobalTrueVSize();
  if (myid == 0) {
    cout << "Number of global H(Div) unknowns: " << size_rt << endl;
  }
  if (myid == 0) {
    cout << "Number of local L2 unknowns: " << L2FESpace.GetTrueVSize() << endl;
  }
  // Create the Electromagnetic solver
  MagnetodynamicSolver MagnetEM(HCurlFESpace, HDivFESpace, L2FESpace, ess_bdr_,
                                ess_domain_, muInv, ItrFunc, H_BCFunc);
  if (myid == 0) {
    cout << "Starting initialization Magnet solver." << endl;
  }

  // Display the current number of DoFs in each finite element space
  MagnetEM.PrintSizes();
  // Set the initial conditions for both the current density J and magnetic
  // fields H
  MagnetEM.SetInitial_Hfield();
  // Set the largest stable time step
  double dtmax = 0.001;

  // Create the ODE solver
  BackwardEulerSolver BESolver;
  BESolver.Init(MagnetEM);
  if (myid == 0) {
    cout << "Initialization ODE solver finished." << endl;
  }
  // Initialize VisIt visualization
  VisItDataCollection visit_dc("Magnetodynamic_newton ", pmesh);
  double t = ti;
  MagnetEM.SetTime(t);

  if (visit) {
    MagnetEM.RegisterVisItFields(visit_dc);
  }
  // Write initial fields to disk for VisIt
  if (visit) {
    // MagnetEM.WriteVisItFields(0);
  }
  // The main time evolution loop.
  int it = 0;
  t = 0;
  while (t < tf) {
    // Run the simulation until a snapshot is needed
    BESolver.Step(H1, t, dt); // Step() includes t += dt     H = H +dHdt*dt

    // Update local DoFs with current true DoFs
    // MagnetEM.SyncGridFuncs();

    // Write fields to disk for VisIt
    if (visit) {
      //  MagnetEM.WriteVisItFields(it);
    }
    it++;
  }

  ParDiscreteCurlOperator *curl =
      new ParDiscreteCurlOperator(&HCurlFESpace, &HDivFESpace);
  curl->Assemble();
  curl->Finalize();
  ParGridFunction H1_(&HCurlFESpace);
  ParGridFunction J1_(&HDivFESpace);
  ParGridFunction J_(&HCurlFESpace);
  ParGridFunction J2_(&L2FESpace);
  H1_.Distribute(H1);
  curl->Mult(H1_, J1_);

  ParMixedBilinearForm *hDivHCurlMuInv_ =
      new ParMixedBilinearForm(&HDivFESpace, &HCurlFESpace);
  hDivHCurlMuInv_->AddDomainIntegrator(new VectorFEMassIntegrator);
  hDivHCurlMuInv_->Assemble();
  hDivHCurlMuInv_->Finalize();
  // VectorGridFunctionCoefficient J_coeff(&J1_);
  // ParDiscreteLinearOperator *Rho_DivCurl_Mass = new
  // ParDiscreteLinearOperator(&HDivFESpace, &L2FESpace);
  // Rho_DivCurl_Mass->AddDomainInterpolator(new
  // VectorInnerProductInterpolator(J_coeff)); Rho_DivCurl_Mass->Assemble();
  // Rho_DivCurl_Mass->Finalize();
  // Rho_DivCurl_Mass->Mult(J1_,J2_ );
  // H1_.ProjectVectorFieldOn(J2,1);
  ostringstream mesh_name, sol_name;
  mesh_name << "meshJ4." << setfill('0') << setw(6) << myid;
  sol_name << "solJ4." << setfill('0') << setw(6) << myid;

  ofstream mesh_ofs(mesh_name.str().c_str());
  mesh_ofs.precision(8);
  pmesh->Print(mesh_ofs);

  ofstream EMnewton_sol_ofs(sol_name.str().c_str());
  EMnewton_sol_ofs.precision(8);
  J1_.Save(EMnewton_sol_ofs);
  // delete Rho_DivCurl_Mass;
  delete hDivHCurlMuInv_;
  delete curl;
  delete fec_L2;
  delete fec_ND;
  delete fec_RT;
  delete pmesh;
  // MFEMFinalizePetsc();
  MPI_Finalize();
  return 0;
}

/** Auxiliary class to provide preconditioners for matrix-free methods */
// class PreconditionerFactory : public PetscPreconditionerFactory
//{
// private:
// const ReducedSystemOperator& op; // unused for now (generates warning)

// public:
// PreconditionerFactory(const Hform_Operator& op_, const string& name_)
//  : PetscPreconditionerFactory(name_) /* , op(op_) */ {}
// virtual mfem::Solver* NewPreconditioner(const mfem::OperatorHandle&);
// virtual ~PreconditionerFactory() {}
//};
// This method gets called every time we need a preconditioner "oh"
// contains the PetscParMatrix that wraps the operator constructed in
// the GetGradient() method (see also PetscSolver::SetJacobianType()).
// In this example, we just return a customizable PetscPreconditioner
// using that matrix. .
// Solver* PreconditionerFactory::NewPreconditioner(const mfem::OperatorHandle&
// oh)
//{
// PetscParMatrix *pP;
// oh.Get(pP);
// cout << " use PCfactory " << endl;
// return new PetscPreconditioner(*pP,"jfnk_"); ////
//}
// the Linearform Integrator to perform the normal flux on the an internal
// interface, used to impose current
class SurfaceFluxLFIntegrator : public LinearFormIntegrator {
  Vector shape;
  VectorCoefficient &Q;
  int oa, ob;

public:
  /// Constructs a boundary integrator with a given Coefficient QG
  SurfaceFluxLFIntegrator(VectorCoefficient &QG, int a = 1, int b = 1)
      : Q(QG), oa(a), ob(b) {}

  virtual void AssembleRHSElementVect(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      Vector &elvect);
  using LinearFormIntegrator::AssembleRHSElementVect;
};
void SurfaceFluxLFIntegrator::AssembleRHSElementVect(const FiniteElement &el,
                                                     ElementTransformation &Tr,
                                                     Vector &elvect) {
  int dim = el.GetDim() + 1;
  int dof = el.GetDof();
  Vector nor(dim), Qvec;
  shape.SetSize(dof);
  elvect.SetSize(dof);
  elvect = 0.0;

  const IntegrationRule *ir = IntRule;
  if (ir == NULL) {
    int intorder = oa * el.GetOrder() + ob; // <----------
    ir = &IntRules.Get(el.GetGeomType(), intorder);
  }

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    Tr.SetIntPoint(&ip);
    CalcOrtho(Tr.Jacobian(), nor);
    Q.Eval(Qvec, Tr, ip);

    el.CalcPhysShape(Tr, shape);
    double area = nor.Norml2();
    double prod = Qvec * nor;
    double signedArea =
        area * ((fabs(prod) < 1e-4 * fabs(area)) ? 0.0 : copysign(1.0, prod));
    elvect.Add(ip.weight * signedArea, shape);
  }
}

class SuperconductorEJIntegrator : public NonlinearFormIntegrator {
protected:
  Coefficient *Q;

private:
  DenseMatrix curlshape, curlshape_dFt, Jac;
  Vector vec, pointflux;
  double E0 = 0.001;
  int n = 20;
  double Jc = 20000000;

public:
  SuperconductorEJIntegrator(Coefficient &m) : Q(&m) {}
  virtual void AssembleElementGrad(const FiniteElement &el,
                                   ElementTransformation &Ttr,
                                   const Vector &elfun, DenseMatrix &elmat);
  virtual void AssembleElementVector(const FiniteElement &el,
                                     ElementTransformation &Ttr,
                                     const Vector &elfun, Vector &elvect);
};

void SuperconductorEJIntegrator::AssembleElementGrad(const FiniteElement &el,
                                                     ElementTransformation &Ttr,
                                                     const Vector &elfun,
                                                     DenseMatrix &elmat) {
  int nd = el.GetDof(), dim = el.GetDim();
  Vector rho(dim), J(dim), curl_dot_trial(nd), curl_dot_test(nd),
      corrector_vec(nd);
  DenseMatrix curlshape(nd, dim), curlshape_test(nd, dim),
      curlshape_dFt(nd, dim), curlphyshape(nd, dim),
      elmat_com(nd); // both trial and test function in Nedelec space,
                     // represented with curlshape
  elmat.SetSize(nd);
  elmat = 0.0;
  double w;
  corrector_vec[0] = 1e-99;
  corrector_vec[1] = 1e-99;
  corrector_vec[2] = 1e-99;
  corrector_vec[3] = 1e-99;
  corrector_vec[4] = 1e-99;
  corrector_vec[5] = 1e-99;
  const IntegrationRule *ir = IntRule;
  if (ir == NULL) {
    int order;
    if (el.Space() == FunctionSpace::Pk) {
      order = 2 * el.GetOrder() - 2;
    } else {
      order = 2 * el.GetOrder();
    }
    ir = &IntRules.Get(el.GetGeomType(), order);
  }
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    Ttr.SetIntPoint(&ip);
    w = ip.weight / Ttr.Weight();
    el.CalcCurlShape(ip, curlshape); // curl operation on the shape function
    el.CalcPhysCurlShape(Ttr, curlphyshape);
    MultABt(
        curlshape, Ttr.Jacobian(),
        curlshape_dFt); // the curl operation of H(curl) basis space, (mapped to
                        // the reference element): H(div)  u(x) = (J/w) * uh(xh)

    // AddMult_a_AAt(1.0, curlshape_dFt, elmat_com);// (Curl u, curl v)*J_de*w
    // cout << " elmat " << elmat_com(1,0) << elmat_com(1,1)  <<elmat_com(1,2)
    // << elmat_com(1,3) <<elmat_com(1,4) << elmat_com(1,5) <<endl;

    if ((Q->Eval(Ttr, ip)) == 0) {
      // double J_norm = J.Norml2();
      // double J_de = 0.1*(E0/Jc); //*n*pow(((J_norm)/Jc), n-1); // derivative
      // factor (n-1)*E0/Jc*(CurlH.Norm/Jc)^(n-2) need to change to n
      // w *= J_de;
      // curl_dot_trial *= n*(E0/Jc)*(1/(1e-99 + pow(J_norm,
      // 2)))*pow(((J_norm)/Jc), n)*w; AddMultVWt(curl_dot_trial, curl_dot_test,
      // elmat);
      //}
      curlphyshape.MultTranspose(elfun, J);
      // double J2 = Q_J->Eval(Ttr, ip);
      double J_norm = J.Norml2();
      double J_de = n * (E0 / Jc) * pow(((J_norm) / Jc), n - 1);
      w *= J_de;
      // w *= 0.00000001;
      // cout << J.Norml2() << "  J norm  " << J_norm << endl;
      AddMult_a_AAt(w, curlshape_dFt, elmat);
      curlshape_dFt.Mult(J, curl_dot_test);
      curlshape_dFt.Mult(J, curl_dot_trial);
      add(curl_dot_trial, corrector_vec, curl_dot_trial);
      double a = (n - 1) * (1 / (1e-99 + pow(J_norm, 2))) * (E0 / Jc) *
                 pow(((J_norm) / Jc), n - 1) * (ip.weight / Ttr.Weight());
      AddMult_a_VWt(a, curl_dot_trial, curl_dot_test, elmat);
    } else {
      w *= Q->Eval(Ttr, ip);
      AddMult_a_AAt(w, curlshape_dFt, elmat);
    }
  }
}

void SuperconductorEJIntegrator::AssembleElementVector(
    const FiniteElement &el, ElementTransformation &Ttr, const Vector &elfun,
    Vector &elvect) {
  int nd = el.GetDof(), dim = el.GetDim();
  DenseMatrix curlshape(nd, dim), curlshape_dFt(nd, dim), elmat_te(nd),
      curlphyshape(nd, dim); // both trial and test function in Nedelec space,
                             // represented with curlshape
  elmat_te = 0.0;
  double w;
  Vector rho(dim), J(dim);
  elvect.SetSize(nd);
  elvect = 0.0;
  // cout << "elfun size " <<  elfun.Size() << endl;
  // cout << "Densemtrix row col " << nd <<" Elfun size " << elfun.Size() <<
  // endl;
  const IntegrationRule *ir = IntRule;
  if (ir == NULL) {
    int order;
    if (el.Space() == FunctionSpace::Pk) {
      order = 2 * el.GetOrder() - 2;
    } else {
      order = 2 * el.GetOrder();
    }
    ir = &IntRules.Get(el.GetGeomType(), order);
  }
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const IntegrationPoint &ip = ir->IntPoint(i);
    Ttr.SetIntPoint(&ip);
    w = ip.weight / Ttr.Weight();
    el.CalcCurlShape(ip, curlshape); // curl operation on the shape function

    el.CalcPhysCurlShape(Ttr, curlphyshape);
    MultABt(curlshape, Ttr.Jacobian(), curlshape_dFt);
    curlphyshape.MultTranspose(
        elfun,
        J); // compute the current density J, elfun is the vector values of H
    // if ((Q->Eval(Ttr, ip)) == 0) {cout << "J value " <<  J_norm << endl; }
    double J_norm = J.Norml2();
    if ((Q->Eval(Ttr, ip)) == 0) {
      // cout << J_norm << " J norm  " << J.Norml2() << endl;
      double J_f = (E0 / Jc) * pow(((J_norm) / Jc), n - 1);
      w *= J_f;
      // w *= 0.00000001;
    } else {
      w *= Q->Eval(Ttr, ip);
    }
    AddMult_a_AAt(w, curlshape_dFt, elmat_te);
    elmat_te.AddMult(elfun, elvect);
  }
}
class Hform_Operator : public Operator {
private:
  ParBilinearForm *LHS_linear, *LHS_massH, *General_mass;
  // ParDiscreteLinearOperator  *Rho_DivCurl_Mass;
  ParNonlinearForm *LHS_nonlinear, *Rho_nonlinear;
  ParFiniteElementSpace *HCurlFESpace, *HDivFESpace, *L2FESpace;
  mutable HypreParMatrix *Jacobian;
  double dt;
  const Vector *H_rhs, *currents;
  const Array<int> ess_tdof_list;
  Coefficient *rhoCoef_L;
  int io;
  int myid_;

public:
  Hform_Operator(MPI_Comm comm, int myid_, ParBilinearForm *LHS_linear_,
                 ParBilinearForm *LHS_massH_, ParBilinearForm *General_mass_,
                 ParFiniteElementSpace *HCurlFESpace_,
                 ParFiniteElementSpace *HDivFESpace_,
                 ParFiniteElementSpace *L2FESpace_,
                 const Array<int> *ess_tdof_list_, Coefficient *rhoCoef_L_);
  /// Set current dt, v, x values - needed to compute action and Jacobian.
  void SetParameters(double dt_, const Vector *H_);
  /// Compute y = H(x + dt (v + dt k)) + M k + S (v + dt k).
  virtual void Mult(const Vector &k, Vector &Residual) const;
  /// Compute J = M + dt S + dt^2 grad_H(x + dt (v + dt k)).
  virtual Operator &GetGradient(const Vector &k) const;
  virtual ~Hform_Operator();
};
//
Hform_Operator::Hform_Operator(
    MPI_Comm comm, int myid_, ParBilinearForm *LHS_linear_,
    ParBilinearForm *LHS_massH_, ParBilinearForm *General_mass_,
    ParFiniteElementSpace *HCurlFESpace_, ParFiniteElementSpace *HDivFESpace_,
    ParFiniteElementSpace *L2FESpace_, const Array<int> *ess_tdof_list_,
    Coefficient *rhoCoef_L_)
    : Operator((myid_ == 0) ? LHS_linear_->ParFESpace()->TrueVSize()
                            : LHS_linear_->ParFESpace()->TrueVSize()),
      LHS_nonlinear(NULL), LHS_linear(LHS_linear_), HCurlFESpace(HCurlFESpace_),
      HDivFESpace(HDivFESpace_), L2FESpace(L2FESpace_), LHS_massH(LHS_massH_),
      General_mass(General_mass_), Rho_nonlinear(NULL), dt(0.0),
      ess_tdof_list(*ess_tdof_list_), rhoCoef_L(rhoCoef_L_), Jacobian(NULL),
      myid_(0) {
  // MPI_Comm_rank(LHS_linear_->ParFESpace()->GetComm(), &myid_);

  LHS_nonlinear = new ParNonlinearForm(
      HCurlFESpace); // H_formulation left hand side, the nonlinear component
  LHS_nonlinear->AddDomainIntegrator(new SuperconductorEJIntegrator(
      *rhoCoef_L)); // the curl_curl operation with nonlinear E-J
  LHS_nonlinear->SetEssentialTrueDofs(ess_tdof_list);

  if (myid_ == 0) {
    cout << "Hform constructed " << endl;
  }
}

void Hform_Operator::SetParameters(double dt_, const Vector *H_) {
  dt = dt_;
  H_rhs = H_;
}

void Hform_Operator::Mult(const Vector &k,
                          Vector &Residual) const // k should be the trueDofs
{
  Vector w(k.Size());
  Residual.SetSize(k.Size());
  // cout << "Mult start " << endl;
  LHS_nonlinear->Mult(k, Residual);
  // cout << "nonlinear norm " << Residual.Norml2() << endl;
  // Residual =0.0;
  // LHS_linear->TrueAddMult(k, Residual);
  add(k, -1, *H_rhs, w);
  w *= 1.0 / dt;
  LHS_massH->TrueAddMult(w, Residual);
  if (myid_ == 0) {
    cout << "nonlinear Mult updated with mass" << endl;
  }

  Residual.SetSubVector(ess_tdof_list, 0.0);
  // cout << "total  norm " << Residual.Norml2() << endl;
  if (myid_ == 0) {
    cout << " Mult finished" << endl;
  }
}

Operator &Hform_Operator::GetGradient(const Vector &k) const {
  delete Jacobian;
  SparseMatrix *localJ = Add(0.0, LHS_linear->SpMat(), 1.0 / dt,
                             LHS_massH->SpMat()); // curlcurl(H)+H/dt

  LHS_nonlinear->SetEssentialTrueDofs(ess_tdof_list);
  // cout << "Jacobian start " << endl;
  // SparseMatrix *ABt;
  // ABt->SetSize(k.Size())

  localJ->Add(1.0, LHS_nonlinear->GetLocalGradient(k));
  // ABt = mfem::Mult( LHS_linear->SpMat(), LHS_nonlinear->GetLocalGradient(k),
  // NULL);
  // if we are using PETSc, the HypreParCSR Jacobian will be converted to
  // PETSc's AIJ on the fly
  Jacobian = General_mass->ParallelAssemble(localJ);
  delete localJ;
  HypreParMatrix *Je = Jacobian->EliminateRowsCols(ess_tdof_list);
  delete Je;
  if (myid_ == 0) {
    cout << "Jacobian finished " << endl;
  }
  return *Jacobian;
}

Hform_Operator::~Hform_Operator() {
  delete Jacobian;
  delete LHS_nonlinear;
}
MagnetodynamicSolver::MagnetodynamicSolver(
    ParFiniteElementSpace &HCurlFESpace_, ParFiniteElementSpace &HDivFESpace_,
    ParFiniteElementSpace &L2FESpace_, Array<int> &ess_bdr_,
    Array<int> &ess_domain_, double (*muInv)(const Vector &),
    double (*ItrFunc)(double),
    void (*H_BCFunc)(const Vector &, double, Vector &))
    : HCurlFESpace(&HCurlFESpace_), HDivFESpace(&HDivFESpace_),
      L2FESpace(&L2FESpace_), ess_bdr(ess_bdr_), ess_domain(ess_domain_),
      myid(0), visit_dc_(NULL), muInv_(muInv), H_BC(H_BCFunc),
      ItrFunc_(ItrFunc), H_t(NULL), H_BCCoef_(NULL), H_t1(NULL), J_t1(NULL),
      curl_(NULL), Hform_LHS_linear(NULL), Hform_LHS_nonlinear(NULL),
      Hform_massH(NULL), General_mass(NULL), ess_tdof_list(NULL),
      Hform_oper(NULL), newton_solver(HCurlFESpace->GetComm())

{
  MPI_Comm_rank(HCurlFESpace->GetComm(), &myid);
  if (myid == 0) {
    cout << "Creating Magnet Solver" << endl;
  }
  const double rel_tol = 0.02;
  // Select surface attributes for Dirichlet BCs
  Array<int> ess_tdof_list_, ess_bdr_cs(ess_bdr.Size());
  int brs = ess_bdr.Size();
  int dms = ess_domain.Size();
  if (myid == 0) {
    cout << "Number of domains " << dms << endl;
  }
  if (myid == 0) {
    cout << "Number of boundaries " << brs << endl;
  }
  ess_bdr = 1.0;
  // ess_bdr[ess_bdr.Size()-1] = 0;// the last surface attribute is the internal
  // interface (cross section)
  // ess_bdr_ makes the outer boundary, ess_bdr_cs marks the inner interface
  // (cross section)
  ess_bdr_cs = 0;
  ess_bdr_cs[ess_bdr.Size() - 1] = 1;
  // Setup various coefficients
  HCurlFESpace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list_);
  ess_tdof_list = &ess_tdof_list_;
  if (myid == 0) {
    cout << "essential true dofs num " << ess_tdof_list->Size() << endl;
  }
  // magnetic field H on the outer surface   define a coefficient from function
  // define the time-dependent coefficient for the BC, with the BC function
  Vector rho_L(ess_domain.Size()); // resistivity for the linear domains,
                                   // superconducting domains set to 0
  rho_L(1) = 0.0;
  rho_L(0) = 1.0 * 100;
  rhoCoef_L = new PWConstCoefficient(rho_L);
  Vector rho_NL(ess_domain.Size()); // resistivity for nonlinear domains,
                                    // non-superconducting domains set to 0
  // rho_NL= 0.0;
  rho_NL(1) = 1.0;
  rho_NL(0) = 0.0;
  rhoCoef_NL = new PWConstCoefficient(rho_NL);

  muCoef_ = new FunctionCoefficient(muInv_); // coefficient for permeability
  H_BCCoef_ = new VectorFunctionCoefficient(
      3, H_BC); // coefficient for time-dependent boundary condition for H

  if (myid == 0) {
    cout << "Creating Coefficients" << endl;
  }
  const int skip_zero_entries = 0;
  const double ref_density = 1.0; // density in the reference configuration
  ConstantCoefficient ones_coef(ref_density);
  // Bilinear Forms
  Hform_LHS_linear = new ParBilinearForm(
      HCurlFESpace); // H_formulation left hand side, the linear component
  Hform_LHS_linear->AddDomainIntegrator(
      new CurlCurlIntegrator(*rhoCoef_L)); // the curl_curl operation on LHS of
                                           // the linear part curl(curl H2)
  Hform_LHS_linear->Assemble(skip_zero_entries);
  Hform_LHS_linear->Finalize(skip_zero_entries);
  if (myid == 0) {
    cout << "Creating bilinear objects" << endl;
  }

  Hform_massH = new ParBilinearForm(
      HCurlFESpace); ////The mass operation, compute the H2/dt
  Hform_massH->AddDomainIntegrator(new VectorFEMassIntegrator(*muCoef_));
  Hform_massH->Assemble(skip_zero_entries);
  Hform_massH->Finalize(skip_zero_entries);

  // Hform_LHS_nonlinear = new ParNonlinearForm (HCurlFESpace); // H_formulation
  // left hand side, the nonlinear component
  // Hform_LHS_nonlinear->AddDomainIntegrator(new
  // SuperconductorEJIntegrator(*rhoCoef_L, *rhoCoef_NL )); // the curl_curl
  // operation with nonlinear E-J
  // Hform_LHS_nonlinear->SetEssentialTrueDofs(*ess_tdof_list);
  General_mass = new ParBilinearForm(HCurlFESpace);
  General_mass->AddDomainIntegrator(new VectorFEMassIntegrator(ones_coef));
  General_mass->Assemble(skip_zero_entries);
  General_mass->Finalize(skip_zero_entries);

  curl_ = new ParDiscreteCurlOperator(HCurlFESpace, HDivFESpace);
  curl_->Assemble();
  curl_->Finalize();
  H_t1 = new ParGridFunction(HCurlFESpace); // the solution vector for H field
  J_t1 = new ParGridFunction(HDivFESpace);  // Current density

  //.......................................................................................
  // for total net transport currents in each conductor
  // use PETS nonlinear solver
  Hform_oper = new Hform_Operator(
      MPI_COMM_WORLD, myid, Hform_LHS_linear, Hform_massH, General_mass,
      HCurlFESpace, HDivFESpace, L2FESpace, ess_tdof_list, rhoCoef_L);
  if (myid == 0) {
    cout << "H form Operator created" << endl;
  }
  cout << "comm is " << HCurlFESpace->GetComm() << endl;

  HypreSmoother *J_hypreSmoother = new HypreSmoother;
  J_hypreSmoother->SetType(HypreSmoother::l1Jacobi);
  J_hypreSmoother->SetPositiveDiagonal(true);
  J_prec = J_hypreSmoother;

  MINRESSolver *J_minres = new MINRESSolver(HCurlFESpace->GetComm());
  J_minres->SetRelTol(rel_tol);
  J_minres->SetAbsTol(0.0);
  J_minres->SetMaxIter(300);
  J_minres->SetPrintLevel(-1);
  J_minres->SetPreconditioner(*J_prec);
  J_solver = J_minres;

  newton_solver.iterative_mode = true;
  newton_solver.SetSolver(*J_solver);
  newton_solver.SetOperator(*Hform_oper);
  newton_solver.SetPrintLevel(1); // print Newton iterations
  newton_solver.SetRelTol(rel_tol);
  newton_solver.SetAbsTol(0.0);
  newton_solver.SetAdaptiveLinRtol(2, 0.5, 0.9);
  newton_solver.SetMaxIter(300);
  //............................................................................................
}

MagnetodynamicSolver::~MagnetodynamicSolver() {
  delete J_solver;
  delete J_prec;

  delete Hform_oper;
  delete J_t1;
  delete H_t1;
  // delete current_constraint;
  delete muCoef_;
  delete H_BCCoef_;
  delete rhoCoef_NL;
  delete rhoCoef_L;
  delete curl_;
  delete Hform_LHS_linear;
  delete Hform_LHS_nonlinear;
  delete General_mass;
  delete Hform_massH;
  // delete H_t;
}

HYPRE_Int MagnetodynamicSolver::GetProblemSize() {
  return HCurlFESpace->GlobalTrueVSize();
}

void MagnetodynamicSolver::PrintSizes() {
  HYPRE_Int size_nd = HCurlFESpace->GlobalTrueVSize();
  HYPRE_Int size_rt = HDivFESpace->GlobalTrueVSize();
  cout << "Number of H(Curl) unknowns: " << size_nd << endl;
  cout << "Number of H(Div)  unknowns: " << size_rt << endl;
}

void MagnetodynamicSolver::SetInitial_Hfield() {
  *H_t1 = 0.0;
  H_t = H_t1->ParallelProject(); // initialize *H_t=0
  cout << "H field initialized " << endl;
}

void MagnetodynamicSolver::ImplicitSolve(
    double dt, const Vector &H,
    Vector &dHdt) // perform the calculation of H and J in every time step
{
  if (myid == 0) {
    cout << "Time t= " << t << endl;
  }
  dHdt.SetSize(H.Size());
  dHdt = 0.0;
  H_BCCoef_->SetTime(t); // update the BC of H to current time step, 't' is
                         // member data from mfem::TimeDependentOperator
  H_t2 = new ParGridFunction(HCurlFESpace); // the solution vector for H field
  *H_t2 = 0.0; // initialize the solution vector for the current time step to 0
  H_t2->ProjectBdrCoefficientTangent(
      *H_BCCoef_, ess_bdr); //    add Dirichlet BC into solution vector
  Vector H_rhs(H.Size());   // empty vector
  H_rhs = H;
  HypreParVector *H_par =
      H_t2->GetTrueDofs(); // initilized the solution vector with BC
  Vector H_lhs(H.Size());
  H_lhs = *H_par;

  if (myid == 0) {
    cout << "BC values initialized " << H_par->Norml2() << endl;
  }

  // Vector h_lhs= *H_par->GlobalVector();; // the x vector LHS

  Hform_oper->SetParameters(dt, &H_rhs);
  Vector zero;
  newton_solver.Mult(zero, H_lhs); // Mult(b,x);
  MFEM_VERIFY(newton_solver.GetConverged(), "Newton solver did not converge.");
  H_t1->Distribute(H_lhs);
  curl_->Mult(*H_t1, *J_t1);
  // cout << "current density " << J_t1->Norml2() << endl;
  add(H_lhs, -1.0, H, dHdt);
  dHdt *= 1 / dt;
  if (myid == 0) {
    cout << "time step finished " << endl;
  }
}
void MagnetodynamicSolver::SyncGridFuncs() {
  H_t1->Distribute(*H_t); // H_ returned by ODE_solver->step()
}

void MagnetodynamicSolver::RegisterVisItFields(VisItDataCollection &visit_dc) {
  visit_dc_ = &visit_dc;

  // visit_dc.RegisterField("H", H_t1);
  visit_dc.RegisterField("J", J_t1);
  // visit_dc.RegisterField("B", B_);
}

void MagnetodynamicSolver::WriteVisItFields(int it) {
  if (visit_dc_) {
    if (myid == 0) {
      cout << "Writing VisIt files ..." << flush;
    }

    HYPRE_Int prob_size = this->GetProblemSize();
    visit_dc_->SetCycle(it);
    visit_dc_->SetTime(prob_size);
    visit_dc_->Save();

    if (myid == 0) {
      cout << " done." << endl;
    }
  }
}
