// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "Joint_solver.hpp"

#ifdef MFEM_USE_MPI

using namespace std;

namespace mfem
{

using namespace common;

namespace electromagnetics
{

double PWGridfunctionCoefficient::Eval(ElementTransformation &T,
                                      const IntegrationPoint &ip)
{
   // given the function and the GridFunction, produce a nonlinear GridFunctionCoefficient
   Mesh *gf_mesh = GridF->FESpace()->GetMesh();
   double field = GridF->GetValue(T, ip, Component); 
   return Function(field) ;
}

Joint_solver::Joint_solver (
   int stateVectorLen,
   ParFiniteElementSpace &L2FES,
   ParFiniteElementSpace &HCurlFES,
   ParFiniteElementSpace &HDivFES,
   ParFiniteElementSpace &HGradFES,
   Array<int> &domains_arg,  
   Array<int> &ess_bdr_arg,
   Array<int> &thermal_ess_bdr_arg,                 
   Array<int> &poisson_ess_bdr_arg, Array<int> &poisson_ess_bdr_in_arg, Array<int> &poisson_ess_bdr_out_arg, Array<int> &poisson_ess_bdr_outer_arg, 
   Array<int> &magnetic_nat_bdr_x_arg, Array<int> &magnetic_nat_bdr_y_arg, Array<int> &magnetic_nat_bdr_z_arg,
   double mu_coef,
   std::map<int, double> sigmaAttMap,
   std::map<int, double> TcapacityAttMap,
   std::map<int, double> InvTcapAttMap,
   std::map<int, double> InvTcondAttMap,
   double  (*ItrFunc_ ) (double),
   void   (*dBdt_xyz_ )(const Vector&, double, Vector&),
   double  (*T_sigma_copper_ ) (double), double (*T_sigma_HTS_ ) (double),double  (*T_sigma_X_) (double)                      
   )

   : TimeDependentOperator(stateVectorLen, 0.0),
     L2FESpace(L2FES), HCurlFESpace(HCurlFES), HDivFESpace(HDivFES),
     HGradFESpace(HGradFES), myid(0),
     a0(NULL), a1(NULL), a2(NULL), m1(NULL), m2(NULL), m3(NULL),
     s1(NULL), s2(NULL), grad(NULL), curl(NULL), weakDiv(NULL), weakDivC(NULL),
     weakCurl(NULL),Ecurl_Jdiv_(NULL),Jdiv_Jdiv_(NULL), Ecurl_Jcurl_(NULL),
     A0(NULL), A1(NULL), A2(NULL), M1(NULL), M2(NULL), M3(NULL),
     X0(NULL), X1(NULL), X2(NULL), B0(NULL), B1(NULL), B2(NULL), B3(NULL),
     v0(NULL), v1(NULL), v2(NULL),Jcurl_(NULL),Phi_rhs_(NULL),
     amg_a0(NULL), pcg_a0(NULL), ads_a2(NULL), pcg_a2(NULL), ams_a1(NULL),
     pcg_a1(NULL), dsp_m3(NULL),pcg_m3(NULL),
     dsp_m1(NULL), pcg_m1(NULL), dsp_m2(NULL), pcg_m2(NULL),
     mu(mu_coef), dt_A1(-1.0), dt_A2(-1.0),
     ItrFunc(ItrFunc_),dBdt_xyz(dBdt_xyz_), dBdt_Coeff_xyz(NULL),
     T_sigma_copper(T_sigma_copper_), T_sigma_HTS(T_sigma_HTS_), T_sigma_X(T_sigma_X_)

{
   MPI_Comm_rank(HDivFESpace.GetComm(), &myid);
   ess_bdr.SetSize(ess_bdr_arg.Size());
   domains.SetSize(domains_arg.Size());
   for (int i=0; i<ess_bdr_arg.Size(); i++)
   {
      ess_bdr[i] = ess_bdr_arg[i];
   }
   thermal_ess_bdr.SetSize(thermal_ess_bdr_arg.Size());
   for (int i=0; i<thermal_ess_bdr_arg.Size(); i++)
   {
      thermal_ess_bdr[i] = thermal_ess_bdr_arg[i];
   }
   poisson_ess_bdr.SetSize(poisson_ess_bdr_arg.Size());
   for (int i=0; i<poisson_ess_bdr_arg.Size(); i++)
   {
      poisson_ess_bdr[i] = poisson_ess_bdr_arg[i];
   }
   poisson_ess_bdr_in.SetSize(poisson_ess_bdr_in_arg.Size());
   for (int i=0; i<poisson_ess_bdr_in_arg.Size(); i++)
   {
      poisson_ess_bdr_in[i] = poisson_ess_bdr_in_arg[i];
   }
   poisson_ess_bdr_out.SetSize(poisson_ess_bdr_out_arg.Size());
   for (int i=0; i<poisson_ess_bdr_out_arg.Size(); i++)
   {
      poisson_ess_bdr_out[i] = poisson_ess_bdr_out_arg[i];
   }
   poisson_ess_bdr_outer.SetSize(poisson_ess_bdr_outer_arg.Size());
   for (int i=0; i<poisson_ess_bdr_outer_arg.Size(); i++)
   {
      poisson_ess_bdr_outer[i] = poisson_ess_bdr_outer_arg[i];
   }
   magnetic_nat_bdr_x.SetSize(magnetic_nat_bdr_x_arg.Size());
   for (int i=0; i<magnetic_nat_bdr_x_arg.Size(); i++)
   {
      magnetic_nat_bdr_x[i] = magnetic_nat_bdr_x_arg[i];
   }
   magnetic_nat_bdr_y.SetSize(magnetic_nat_bdr_y_arg.Size());
   for (int i=0; i<magnetic_nat_bdr_y_arg.Size(); i++)
   {
      magnetic_nat_bdr_y[i] = magnetic_nat_bdr_y_arg[i];
   }
   magnetic_nat_bdr_z.SetSize(magnetic_nat_bdr_z_arg.Size());
   for (int i=0; i<magnetic_nat_bdr_z_arg.Size(); i++)
   {
      magnetic_nat_bdr_z[i] = magnetic_nat_bdr_z_arg[i];
   }
   sigma     = new MeshDependentCoefficient(sigmaAttMap);
   Tcapacity = new MeshDependentCoefficient(TcapacityAttMap);
   InvTcap   = new MeshDependentCoefficient(InvTcapAttMap);
   InvTcond  = new MeshDependentCoefficient(InvTcondAttMap);

   
   this->buildM3(*Tcapacity);  
   this->buildM2(*InvTcond);
   this->buildS2(*InvTcap);
   this->buildS1(1.0/mu);
   this->buildCurl(1.0/mu);
   this->buildDiv(*InvTcap);
   this->buildGrad();
     
   dBdt_Coeff_xyz = new VectorFunctionCoefficient(3, *dBdt_xyz); 
   //B_Coeff_z = new VectorFunctionCoefficient(3, *Bt_z); 
   //B_Coeff_x_last = new VectorFunctionCoefficient(3, *Bt_x); 
   //B_Coeff_y_last = new VectorFunctionCoefficient(3, *Bt_y); 
   //B_Coeff_z_last = new VectorFunctionCoefficient(3, *Bt_z); 

   v0 = new ParGridFunction(&HGradFESpace);
   v1 = new ParGridFunction(&HCurlFESpace);
   v2 = new ParGridFunction(&HDivFESpace);
   A0 = new HypreParMatrix;
   A1 = new HypreParMatrix;
   A2 = new HypreParMatrix;
   X0 = new Vector;
   X1 = new Vector;
   X2 = new Vector;
   B0 = new Vector;
   B1 = new Vector;
   B2 = new Vector;
   B3 = new Vector;
  if (myid ==0){ cout << " Joint solver constructed "  << endl; } 
}

void Joint_solver::Init( Vector &X)
{
   Vector zero_vec(3); zero_vec = 0.0;
   VectorConstantCoefficient Zero_vec(zero_vec);
   ConstantCoefficient Zero(0.0);
   ConstantCoefficient T_op(20.0);
    if (myid ==0){ cout << " implicit solve start 11"  << endl; }
   // The big BlockVector stores the fields as follows:
   //    Temperature
   //    Temperature Flux
   //    P field
   //    E field
   //    B field
   //    Joule Heating

   int Vsize_l2 = L2FESpace.GetVSize();
   int Vsize_nd = HCurlFESpace.GetVSize();
   int Vsize_rt = HDivFESpace.GetVSize();
   int Vsize_h1 = HGradFESpace.GetVSize();

   Array<int> true_offset(9);
   true_offset[0] = 0;
   true_offset[1] = true_offset[0] + Vsize_h1;
   true_offset[2] = true_offset[1] + Vsize_l2;
   true_offset[3] = true_offset[2] + Vsize_rt;
   true_offset[4] = true_offset[3] + Vsize_rt;
   true_offset[5] = true_offset[4] + Vsize_nd;
   true_offset[6] = true_offset[5] + Vsize_rt;
   true_offset[7] = true_offset[6] + Vsize_l2;
   true_offset[8] = true_offset[7] + Vsize_rt;

   Vector* xptr = (Vector*) &X;
   ParGridFunction P, E, B, T, F, W, J, J_in;
   P.MakeRef(&HGradFESpace,*xptr,true_offset[0]);
   T.MakeRef(&L2FESpace,   *xptr,true_offset[1]);
   F.MakeRef(&HDivFESpace, *xptr,true_offset[2]);
   J.MakeRef(&HDivFESpace, *xptr,true_offset[3]);
   E.MakeRef(&HCurlFESpace,*xptr,true_offset[4]);
   B.MakeRef(&HDivFESpace, *xptr,true_offset[5]);
   W.MakeRef(&L2FESpace,   *xptr,true_offset[6]);
   J_in.MakeRef(&HDivFESpace,*xptr, true_offset[7]);

   E.ProjectCoefficient(Zero_vec);
   B.ProjectCoefficient(Zero_vec);
   F.ProjectCoefficient(Zero_vec);  
   J.ProjectCoefficient(Zero_vec);
   J_in.ProjectCoefficient(Zero_vec);
   T.ProjectCoefficient(T_op);
   W.ProjectCoefficient(Zero);
   P.ProjectCoefficient(Zero);
   if (myid ==0){ cout << " Vectors Initilisation finished "  << endl; }
}


/*
This is the main computational code that computes dX/dt implicitly
where X is the state vector containing P, E, B, F, T, and W

        S0 P = 0
(M1+dt S1) E = WeakCurl^T B + Grad P
          dB = -Curl E
(M2+dt S2) F = WeakDiv^T T
       M3 dT = WeakDiv F + W
1) solve transport current 
2) solve eddy current 
3) solve joule heating from both 
4) solve temperature
where W is the Joule heating.

Boundary conditions are applied to E.  Boundary conditions are applied to F.  No
boundary conditions are applied to B or T.

The W term in the left hand side is the Joule heating which is a nonlinear
(quadratic) function of E.

P is solution of Div sigma Grad dP = 0.

The total E-field is given by E_tot = E_ind - Grad P, the big equation for E
above is really for E_ind (the induced, or solenoidal, component) and this is
corrected for.
*/
void Joint_solver::ImplicitSolve(const double dt,
                                               const Vector &X, Vector &dX_dt) 
{
   
   if (myid ==0){ cout << " Current time step t = "  << t <<endl; }   
   
   // The big BlockVector stores the fields as follows:
   //    Temperature
   //    Temperature Flux
   //    P field
   //    E field
   //    B field
   //    Joule Heating
  
   int Vsize_l2 = L2FESpace.GetVSize();
   int Vsize_nd = HCurlFESpace.GetVSize();
   int Vsize_rt = HDivFESpace.GetVSize();
   int Vsize_h1 = HGradFESpace.GetVSize();

   Array<int> true_offset(9);
   true_offset[0] = 0;
   true_offset[1] = true_offset[0] + Vsize_h1;
   true_offset[2] = true_offset[1] + Vsize_l2;
   true_offset[3] = true_offset[2] + Vsize_rt;
   true_offset[4] = true_offset[3] + Vsize_rt;
   true_offset[5] = true_offset[4] + Vsize_nd;
   true_offset[6] = true_offset[5] + Vsize_rt;
   true_offset[7] = true_offset[6] + Vsize_l2;
   true_offset[8] = true_offset[7] + Vsize_rt;
   Vector* xptr  = (Vector*) &X;
   
   ParGridFunction P, E, B, T, F, W, J, J_in;
  
   P.MakeRef(&HGradFESpace,*xptr,true_offset[0]);
   T.MakeRef(&L2FESpace,   *xptr,true_offset[1]);
   F.MakeRef(&HDivFESpace, *xptr,true_offset[2]);
   J.MakeRef(&HDivFESpace, *xptr,true_offset[3]);
   E.MakeRef(&HCurlFESpace,*xptr,true_offset[4]);
   B.MakeRef(&HDivFESpace, *xptr,true_offset[5]);
   W.MakeRef(&L2FESpace,   *xptr,true_offset[6]);
   J_in.MakeRef(&HDivFESpace,*xptr, true_offset[7]);

   ParGridFunction dP, dE, dB, dT, dF, dW, dJ,dJ_in;  // here dB means dB/dt
   dP.MakeRef(&HGradFESpace,dX_dt,true_offset[0]);
   dT.MakeRef(&L2FESpace,   dX_dt,true_offset[1]);
   dF.MakeRef(&HDivFESpace, dX_dt,true_offset[2]);
   dJ.MakeRef(&HDivFESpace, dX_dt,true_offset[3]);
   dE.MakeRef(&HCurlFESpace,dX_dt,true_offset[4]);
   dB.MakeRef(&HDivFESpace, dX_dt,true_offset[5]);
   dW.MakeRef(&L2FESpace,   dX_dt,true_offset[6]);
   dJ_in.MakeRef(&HDivFESpace,dX_dt, true_offset[7]);

   dX_dt = 0.0;
     
   PWGridfunctionCoefficient T_sigma1(*T_sigma_X, &T );
   PWGridfunctionCoefficient T_sigma2(*T_sigma_HTS, &T );  
   PWGridfunctionCoefficient T_sigma3(*T_sigma_HTS, &T );
   PWGridfunctionCoefficient T_sigma4(*T_sigma_copper, &T );
   PWGridfunctionCoefficient T_sigma5(*T_sigma_copper, &T );
   Array<Coefficient*> T_sigmaCoeff(domains.Size()); 
   T_sigmaCoeff[0]= &T_sigma1; T_sigmaCoeff[1]= &T_sigma2;  T_sigmaCoeff[2]= &T_sigma3; T_sigmaCoeff[3]= &T_sigma4; T_sigmaCoeff[4]= &T_sigma5;   
   DomainWiseCoefficient T_sigma_Domains(T_sigmaCoeff);

   this->buildA2(*InvTcond, *InvTcap, dt);
   this->buildA1(1.0/mu, T_sigma_Domains, dt);
   this->buildA0(T_sigma_Domains);
   this->buildM1(T_sigma_Domains);
   this->buildA_J(T_sigma_Domains);
   if (myid ==0){ cout << " zero voltage solved "  << endl; }
   // form the Laplacian and solve it
   ParGridFunction Phi_gf(&HGradFESpace); // Electric potential
   Phi_gf = 0.0;
   Array<int> ess_bdr_(poisson_ess_bdr_in.Size()), ess_bdr_in(poisson_ess_bdr_in.Size()), ess_bdr_out(poisson_ess_bdr_in.Size());

   ess_bdr_ = 0; ess_bdr_[4] = 1;ess_bdr_[5] = 1; // both ends of the joint, all essential bdr
   ess_bdr_in = 0; ess_bdr_in[1] = 1; ess_bdr_in[2] = 0;// in for surface 2 (x+)
   ess_bdr_out = 0; ess_bdr_out[1] = 0; ess_bdr_out[2] = 1;// out for surface 3 (z+)
   // the function below is currently not fully supported on AMR meshes
   // Phi_gf.ProjectBdrCoefficient(voltage,poisson_ess_bdr);
   ConstantCoefficient VCoef_in(1.0);
   ConstantCoefficient VCoef_out(0.0);
   ConstantCoefficient VCoef_outer(0.0);
   Phi_gf.ProjectBdrCoefficient(VCoef_in, poisson_ess_bdr_in); 
   Phi_gf.ProjectBdrCoefficient(VCoef_out, poisson_ess_bdr_out); 
   //Phi_gf.ProjectBdrCoefficient(VCoef_outer, poisson_ess_bdr_outer);
   // apply essential BC's and apply static condensation, the new system to
   // solve is A0 X0 = B0
   Array<int> poisson_ess_tdof_list;
   HGradFESpace.GetEssentialTrueDofs(ess_bdr_, poisson_ess_tdof_list);

   *v0 = 0.0;

   a0->FormLinearSystem(poisson_ess_tdof_list,Phi_gf,*Phi_rhs_,*A0,*X0,*B0);

   amg_a0 = new HypreBoomerAMG(*A0); 
  
      pcg_a0 = new HyprePCG(*A0);
      pcg_a0->SetTol(SOLVER_TOL);
      pcg_a0->SetMaxIter(SOLVER_MAX_IT);
      pcg_a0->SetPrintLevel(SOLVER_PRINT_LEVEL);
      pcg_a0->SetPreconditioner(*amg_a0);

   // pcg "Mult" operation is a solve
   // X0 = A0^-1 * B0
   pcg_a0->Mult(*B0, *X0);

   // "undo" the static condensation saving result in grid function dP
   a0->RecoverFEMSolution(*X0,*v0,P);

   if (myid ==0){ cout << " dummy voltage solved "  << endl; } 
   grad->Mult(P,E);
   E *=-1.0; // electric field from the dummy potential for transport currents
   
   if (myid ==0){ cout << " dymmy E field  solved "  << endl; }
   ParGridFunction J_dual(&HDivFESpace);// dual vector of J +
   //E_ = 0.0;
   Ecurl_Jdiv_->Mult(E,J_dual );// compute the dual of dummy J 
   //if ( myid == 0 )  {cout << "calculating the currents" << endl;}
   if (myid ==0){ cout << " dymmy J dual   solved "  << endl; }
   HypreParMatrix MassJdiv;
   Vector j_, jd_;
   ParGridFunction J_(&HDivFESpace);
   Array<int> dbc_dofs_J;
   Jdiv_Jdiv_->FormLinearSystem(dbc_dofs_J, J_, J_dual, MassJdiv, j_, jd_);

   HyprePCG pcgM(MassJdiv);
   pcgM.SetTol(1e-12);
   pcgM.SetMaxIter(500);
   pcgM.SetPrintLevel(0);
   HypreDiagScale diagM;
   pcgM.SetPreconditioner(diagM);
   pcgM.Mult(jd_, j_);
   Jdiv_Jdiv_->RecoverFEMSolution(j_, J_dual, J);// compute the current density J
   if (myid ==0){ cout << " dummy currents solved "  << endl; }
   ParLinearForm I_surface(&HDivFESpace);
   I_surface.AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator, poisson_ess_bdr_in);
   I_surface.Assemble();
   double I_total = I_surface(J);
   if (myid ==0){ cout << " dummy total current is   " << I_total << endl; }
    if (myid ==0){ cout << " real total current is   " << ItrFunc(t) << endl; }
   J *= ItrFunc(t)/I_total;
   E *= ItrFunc(t)/I_total;
   dBdt_Coeff_xyz->SetTime(t); 
   
   this->buildb_B( thermal_ess_bdr, *dBdt_Coeff_xyz, dt);
    // RHS with Neunmann BC ->normal B flux 
   if (myid ==0){ cout << " dBdt at BC created  test all"  << endl; }
   //Ecurl_Jcurl_->AddMult(E, *Jcurl_, -1.0);
   
   ParGridFunction E_dual(&HCurlFESpace);
   E_dual = 0.0;
   //Ecurl_Jcurl_->Mult(E, E_dual);
   m1->Mult(E, E_dual);
   E_dual *=-1.0;
   if (myid ==0){ cout << " transport current subtracted  "  << endl; }
  
   // B is a grid function but weakCurl is not parallel assembled so is OK
   weakCurl->AddMultTranspose(B,E_dual, 1.0); //curlB -sigma E with BC on the RHS (Neumann BC) 
   if (myid ==0){ cout << " E field RHS created  "  << E_dual.Norml2() <<endl; }
   //add(1.0,E_dual,1.0, *Jcurl_,E_dual  );   // use E as a temporary, E = Grad P
   
   if (myid ==0){ cout << " total E rhs created (with BC)"  << E_dual.Norml2() <<endl; }
   ParGridFunction E_in(&HCurlFESpace);
   // form the linear system, including eliminating essential BC's and applying
   // static condensation. The system to solve is A1 X1 = B1
   Array<int> E_ess_tdof_list;
   HCurlFESpace.GetEssentialTrueDofs(poisson_ess_bdr_outer, E_ess_tdof_list);
   E_in = 0.0;
   Array<int> dbc_dofs_E;  // no essential BC for E 
   a1->FormLinearSystem(E_ess_tdof_list,E_in,E_dual,*A1,*X1,*B1);
   if (myid ==0){ cout << " linear system created for E"  << E_dual.Norml2()  <<endl; }
   // We only need to create the solver and preconditioner once
   
      ParFiniteElementSpace *prec_fespace =
         (a1->StaticCondensationIsEnabled() ? a1->SCParFESpace() : &HCurlFESpace);  
      ams_a1 = new HypreAMS(*A1, prec_fespace);
   
      pcg_a1 = new HyprePCG(*A1);
      pcg_a1->SetTol(SOLVER_TOL);
      pcg_a1->SetMaxIter(SOLVER_MAX_IT);
      pcg_a1->SetPrintLevel(SOLVER_PRINT_LEVEL);
      pcg_a1->SetPreconditioner(*ams_a1);
   
   if (myid ==0){ cout << " B1 and X1"   << B1->Norml2()  << X1->Norml2() <<endl; }  
   // solve the system
   // dE = (A1)^-1 [-S1 E]
   pcg_a1->Mult(*B1, *X1);
   if (myid ==0){ cout << " if mult worked can be used "   <<endl; }   
   // this is required because of static condensation, E is a grid function
   a1->RecoverFEMSolution(*X1,E_dual,E_in);
   if (myid ==0){ cout << " induced Efield solved "  << endl; }
   // induced currents computation
   Ecurl_Jdiv_->Mult(E_in,J_dual );// compute the dual of dummy J 
   J_ = 0.0;
   Jdiv_Jdiv_->FormLinearSystem(dbc_dofs_J, J_, J_dual, MassJdiv, j_, jd_);

   HyprePCG pcgM_(MassJdiv);
   pcgM_.SetTol(1e-12);
   pcgM_.SetMaxIter(500);
   pcgM_.SetPrintLevel(0);
   HypreDiagScale diagM_;
   pcgM_.SetPreconditioner(diagM_);
   pcgM_.Mult(jd_, j_);
   Jdiv_Jdiv_->RecoverFEMSolution(j_, J_dual, J_in);// compute the current density J
   if (myid ==0){ cout << " induced currents solved "  <<  J_in.Norml2()<<endl; }
   
   dE = 0.0;
   dJ_in = 0.0;
      
   add(E_in, E, E);// add induced E and transport E
   curl->Mult(E, dB);
   dB *= -1.0;
   dJ = 0.0;
   // Compute Energy Deposition
   this->GetJouleHeating(E,W);
   if (myid ==0){ cout << " Joule heating loss power calculated "  << endl; }
   // v2 = Div^T * W, where W is the Joule heating computed above, and
   // Div is the matrix <div u, v>
   weakDivC->MultTranspose(W, *v2);
   *v2 *= dt;

   // v2 = <v, div u> T + (1.0)*v2
   weakDiv->AddMultTranspose(T, *v2, 1.0);
    if (myid ==0){ cout << " heat flux RHS produced"  << endl; }
   // apply the thermal BC
   Vector zero_vec(3); zero_vec = 0.0;
   VectorConstantCoefficient Zero_vec(zero_vec);
   ParGridFunction F_gf(&HDivFESpace);
   F_gf = 0.0;
   F_gf.ProjectBdrCoefficientNormal(Zero_vec,thermal_ess_bdr);
   if (myid ==0){ cout << " heat flux BC assigned "  << endl; }
   // form the linear system, including eliminating essential BC's and applying
   // static condensation. The system to solve is A2 X2 = B2
   Array<int> thermal_ess_tdof_list;
   HDivFESpace.GetEssentialTrueDofs(thermal_ess_bdr, thermal_ess_tdof_list);
   a2->FormLinearSystem(thermal_ess_tdof_list,F_gf,*v2,*A2,*X2,*B2);
   if (myid ==0){ cout << " Heat flux linear system produced "  << endl; }
   // We only need to create the solver and preconditioner once
  
      ParFiniteElementSpace *prec_fespace1 =
         (a2->StaticCondensationIsEnabled() ? a2->SCParFESpace() : &HDivFESpace);
      ads_a2 = new HypreADS(*A2, prec_fespace1);
      pcg_a2 = new HyprePCG(*A2);
      pcg_a2->SetTol(SOLVER_TOL);
      pcg_a2->SetMaxIter(SOLVER_MAX_IT);
      pcg_a2->SetPrintLevel(SOLVER_PRINT_LEVEL);
      pcg_a2->SetPreconditioner(*ads_a2);
   // solve for dF from a2 dF = v2
   // dF = (A2)^-1 [S2*F + rhs]
   pcg_a2->Mult(*B2, *X2);
   
   // this is required because of static condensation
   a2->RecoverFEMSolution(*X2,*v2,F);
   if (myid ==0){ cout << " Heat flux solved  "  << endl; }
   dF = 0.0;

   // c dT = [W - div F]
   //
   // <u,u> dT = <1/c W,u> - <1/c div v,u>
   //
   // where W is Joule heating and F is the flux that we just computed
   //
   // note: if div is a BilinearForm, then W should be converted to a LoadVector
   // compute load vector <1/c W, u> where W is the Joule heating GF

   // create the Coefficient 1/c W
   //ScaledGFCoefficient Wcoeff(&W, *InvTcap);
   GridFunctionCoefficient Wcoeff(&W);

   // compute <W,u>
   ParLinearForm temp_lf(&L2FESpace);
   temp_lf.AddDomainIntegrator(new DomainLFIntegrator(Wcoeff));
   temp_lf.Assemble();

   // lf = lf - div F
   weakDiv->AddMult(F, temp_lf, -1.0);

   // need to perform mass matrix solve to get temperature T
   // <c u, u> Tdot = -<div v, u> F +  <1/c W, u>
   // NOTE: supposedly we can just invert any L2 matrix, could do that here
   // instead of a solve

    dsp_m3 = new HypreDiagScale(*M3); 
  
      pcg_m3 = new HyprePCG(*M3);
      pcg_m3->SetTol(SOLVER_TOL);
      pcg_m3->SetMaxIter(SOLVER_MAX_IT);
      pcg_m3->SetPrintLevel(SOLVER_PRINT_LEVEL);
      pcg_m3->SetPreconditioner(*dsp_m3);
   

   // solve for dT from M3 dT = lf
   // no boundary conditions on this solve
   pcg_m3->Mult(temp_lf, dT);
   if (myid ==0){ cout << " temperature solved  "  << endl; }
}

void Joint_solver::buildA0(DomainWiseCoefficient &Sigma)
{
   if ( a0 != NULL ) { delete a0; }

   // First create and assemble the bilinear form.  For now we assume the mesh
   // isn't moving, the materials are time independent, and dt is constant. So
   // we only need to do this once.

   // ConstantCoefficient Sigma(sigma);
   ConstantCoefficient zeros(0.0);
   a0 = new ParBilinearForm(&HGradFESpace);
   a0->AddDomainIntegrator(new DiffusionIntegrator(Sigma));
   if (STATIC_COND == 1) { a0->EnableStaticCondensation(); }
   a0->Assemble();
   Phi_rhs_ = new ParLinearForm(&HGradFESpace);
   Phi_rhs_->AddDomainIntegrator(new DomainLFIntegrator(zeros));   
   Phi_rhs_->Assemble();

   // Don't finalize or parallel assemble this is done in FormLinearSystem.
}
void Joint_solver::buildA_J(DomainWiseCoefficient &Sigma)

{
   Ecurl_Jdiv_ = new ParMixedBilinearForm(&HCurlFESpace, &HDivFESpace);
   Ecurl_Jdiv_->AddDomainIntegrator(new VectorFEMassIntegrator(Sigma));
   Ecurl_Jdiv_->Assemble();
   //Ecurl_Jdiv_->Finalize();

   Jdiv_Jdiv_ = new ParBilinearForm( &HDivFESpace);
   Jdiv_Jdiv_->AddDomainIntegrator(new VectorFEMassIntegrator);
   Jdiv_Jdiv_->Assemble();
   //Jdiv_Jdiv_->Finalize();

   Ecurl_Jcurl_ = new ParBilinearForm( &HCurlFESpace);
   Ecurl_Jcurl_->AddDomainIntegrator(new VectorFEMassIntegrator(Sigma));
   Ecurl_Jcurl_->Assemble();
   //Ecurl_Jcurl_->Finalize();
   if (myid ==0){ cout << " currents transformation tools produced  "  << endl; }
}
void Joint_solver::buildb_B( Array<int> &magnetic_nat_bdr, VectorFunctionCoefficient &dBdtCoeff, double dt)
{
  if (myid ==0){ cout << " starting building rhs  with nat BC"  << endl; } 
  if ( Jcurl_ != NULL ) { delete Jcurl_; }
  Vector zero(3); zero = 0.0;
  VectorConstantCoefficient zeros(zero);
  Jcurl_ = new ParLinearForm(&HCurlFESpace);
  Jcurl_->AddDomainIntegrator(new VectorFEDomainLFIntegrator(zeros));   
  //Jcurl_->AddBoundaryIntegrator(new VectorFEBoundaryTangentLFIntegrator(dBdtCoeff_), ess_bdr);
  //Jcurl_->AddBoundaryIntegrator(new VectorFEBoundaryTangentLFIntegrator(dBdtCoeff_y), magnetic_nat_bdr_y);
  Jcurl_->AddBoundaryIntegrator(new VectorFEBoundaryTangentLFIntegrator(dBdtCoeff), magnetic_nat_bdr);
  Jcurl_->Assemble();
  if (myid ==0){ cout << " curlE natural BC specified "  << endl; }
}

void Joint_solver::buildA1(double muInv, DomainWiseCoefficient &Sigma, double dt)
{
   if ( a1 != NULL ) { delete a1; }

   // First create and assemble the bilinear form.  For now we assume the mesh
   // isn't moving, the materials are time independent, and dt is constant. So
   // we only need to do this once.
   
   ConstantCoefficient dtMuInv(dt*muInv);//
   if (myid ==0){ cout << " dtMuiv is "  <<dt*muInv<< endl; }
   a1 = new ParBilinearForm(&HCurlFESpace);
   a1->AddDomainIntegrator(new VectorFEMassIntegrator(Sigma));
   a1->AddDomainIntegrator(new CurlCurlIntegrator(dtMuInv));
   if (STATIC_COND == 1) { a1->EnableStaticCondensation(); }
   a1->Assemble();

   // Don't finalize or parallel assemble this is done in FormLinearSystem.

   dt_A1 = dt;
}

void Joint_solver::buildA2(MeshDependentCoefficient &InvTcond, MeshDependentCoefficient &InvTcap, double dt)
{
   if ( a2 != NULL ) { delete a2; }

   InvTcap.SetScaleFactor(dt);
   a2 = new ParBilinearForm(&HDivFESpace);
   a2->AddDomainIntegrator(new VectorFEMassIntegrator(InvTcond));
   a2->AddDomainIntegrator(new DivDivIntegrator(InvTcap));
   if (STATIC_COND == 1) { a2->EnableStaticCondensation(); }
   a2->Assemble();

   // Don't finalize or parallel assemble this is done in FormLinearSystem.

   dt_A2 = dt;
}

void Joint_solver::buildM1(DomainWiseCoefficient &Sigma)
{
   if ( m1 != NULL ) { delete m1; }

   m1 = new ParBilinearForm(&HCurlFESpace);
   m1->AddDomainIntegrator(new VectorFEMassIntegrator(Sigma));
   m1->Assemble();

   // Don't finalize or parallel assemble this is done in FormLinearSystem.
}

void Joint_solver::buildM2(MeshDependentCoefficient &Alpha)
{
   if ( m2 != NULL ) { delete m2; }

   // ConstantCoefficient MuInv(muInv);
   m2 = new ParBilinearForm(&HDivFESpace);
   m2->AddDomainIntegrator(new VectorFEMassIntegrator(Alpha));
   m2->Assemble();

   // Don't finalize or parallel assemble this is done in FormLinearSystem.
}

void Joint_solver::buildM3(MeshDependentCoefficient &Tcapacity)
{
   if ( m3 != NULL ) { delete m3; }

   // ConstantCoefficient Sigma(sigma);
   m3 = new ParBilinearForm(&L2FESpace);
   m3->AddDomainIntegrator(new MassIntegrator(Tcapacity));
   m3->Assemble();
   m3->Finalize();
   M3 = m3->ParallelAssemble();
}

void Joint_solver::buildS1(double muInv)
{
   if ( s1 != NULL ) { delete s1; }

   ConstantCoefficient MuInv(muInv);
   s1 = new ParBilinearForm(&HCurlFESpace);
   s1->AddDomainIntegrator(new CurlCurlIntegrator(MuInv));
   s1->Assemble();
}

void Joint_solver::buildS2(MeshDependentCoefficient &InvTcap)
{
   if ( s2 != NULL ) { delete s2; }

   // ConstantCoefficient param(a);
   s2 = new ParBilinearForm(&HDivFESpace);
   s2->AddDomainIntegrator(new DivDivIntegrator(InvTcap));
   s2->Assemble();
}

void Joint_solver::buildCurl(double muInv)
{
   if ( curl != NULL ) { delete curl; }
   if ( weakCurl != NULL ) { delete weakCurl; }

   curl = new ParDiscreteLinearOperator(&HCurlFESpace, &HDivFESpace);
   curl->AddDomainInterpolator(new CurlInterpolator);
   curl->Assemble();

   ConstantCoefficient MuInv(muInv);
   weakCurl = new ParMixedBilinearForm(&HCurlFESpace, &HDivFESpace);
   weakCurl->AddDomainIntegrator(new VectorFECurlIntegrator(MuInv));
   weakCurl->Assemble();

   // no ParallelAssemble since this will be applied to GridFunctions
}

void Joint_solver::buildDiv(MeshDependentCoefficient &InvTcap)
{
   if ( weakDiv != NULL ) { delete weakDiv; }
   if ( weakDivC != NULL ) { delete weakDivC; }

   weakDivC = new ParMixedBilinearForm(&HDivFESpace, &L2FESpace);
   weakDivC->AddDomainIntegrator(new VectorFEDivergenceIntegrator(InvTcap));
   weakDivC->Assemble();

   weakDiv = new ParMixedBilinearForm(&HDivFESpace, &L2FESpace);
   weakDiv->AddDomainIntegrator(new VectorFEDivergenceIntegrator());
   weakDiv->Assemble();

   // no ParallelAssemble since this will be applied to GridFunctions
}

void Joint_solver::buildGrad()
{
   if ( grad != NULL ) { delete grad; }

   grad = new ParDiscreteLinearOperator(&HGradFESpace, &HCurlFESpace);
   grad->AddDomainInterpolator(new GradientInterpolator());
   grad->Assemble();

   // no ParallelAssemble since this will be applied to GridFunctions
}

double Joint_solver::ElectricLosses(ParGridFunction &E_gf) const
{
   double el = m1->InnerProduct(E_gf,E_gf);

   double global_el;
   MPI_Allreduce(&el, &global_el, 1, MPI_DOUBLE, MPI_SUM,
                 m2->ParFESpace()->GetComm());

   return el;
}

// E is the input GF, w is the output GF which is assumed to be an L2 scalar
// representing the Joule heating
void Joint_solver::GetJouleHeating(ParGridFunction &E_gf,
                                                 ParGridFunction &w_gf) const
{
   // The w_coeff object stashes a reference to sigma and E, and it has
   // an Eval method that will be used by ProjectCoefficient.
   JouleHeatingCoefficient w_coeff(*sigma, E_gf);

   // This applies the definition of the finite element degrees-of-freedom
   // to convert the function to a set of discrete values
   w_gf.ProjectCoefficient(w_coeff);
}

void Joint_solver::SetTime(const double t_)
{ t = t_; }

Joint_solver::~Joint_solver()
{
   if ( ams_a1 != NULL ) { delete ams_a1; }
   if ( pcg_a1 != NULL ) { delete pcg_a1; }

   if ( dsp_m1 != NULL ) { delete dsp_m1; }
   if ( pcg_m1 != NULL ) { delete pcg_m1; }

   if ( dsp_m2 != NULL ) { delete dsp_m2; }
   if ( pcg_m2 != NULL ) { delete pcg_m2; }

   if ( curl != NULL ) { delete curl; }
   if ( weakDiv != NULL ) { delete weakDiv; }
   if ( weakDivC != NULL ) { delete weakDivC; }
   if ( weakCurl != NULL ) { delete weakCurl; }
   if ( grad != NULL ) { delete grad; }
   if ( Jcurl_ != NULL ) { delete Jcurl_; }
   if ( Ecurl_Jdiv_ != NULL ) { delete Ecurl_Jdiv_; }
   if ( Jdiv_Jdiv_ != NULL ) { delete Jdiv_Jdiv_; }
   if ( Ecurl_Jcurl_ != NULL ) { delete Ecurl_Jcurl_; }
   if ( a0 != NULL ) { delete a0; }
   if ( a1 != NULL ) { delete a1; }
   if ( a2 != NULL ) { delete a2; }
   if ( m1 != NULL ) { delete m1; }
   if ( m2 != NULL ) { delete m2; }
   if ( s1 != NULL ) { delete s1; }
   if ( s2 != NULL ) { delete s2; }

   if ( A0 != NULL ) { delete A0; }
   if ( Phi_rhs_ != NULL ) { delete Phi_rhs_; }  
   if ( X0 != NULL ) { delete X0; }
   if ( B0 != NULL ) { delete B0; }

   if ( A1 != NULL ) { delete A1; }
   if ( X1 != NULL ) { delete X1; }
   if ( B1 != NULL ) { delete B1; }

   if ( A2 != NULL ) { delete A2; }
   if ( X2 != NULL ) { delete X2; }
   if ( B2 != NULL ) { delete B2; }

   if ( v1 != NULL ) { delete v1; }
   if ( v2 != NULL ) { delete v2; }

   if (sigma     != NULL) { delete sigma; }
   if (Tcapacity != NULL) { delete Tcapacity; }
   if (InvTcap   != NULL) { delete InvTcap; }
   if (InvTcond  != NULL) { delete InvTcond; }

   delete amg_a0;
   delete pcg_a0;
   delete pcg_a2;
   delete ads_a2;
   delete m3;
   delete dsp_m3;
   delete pcg_m3;
   delete M1;
   delete M2;
   delete M3;
   delete v0;
   delete B3;
}

double JouleHeatingCoefficient::Eval(ElementTransformation &T,
                                     const IntegrationPoint &ip)
{
   Vector E;
   double thisSigma;
   E_gf.GetVectorValue(T, ip, E);
   thisSigma = sigma.Eval(T, ip);
   return thisSigma*(E*E);
}

MeshDependentCoefficient::MeshDependentCoefficient(
   const std::map<int, double> &inputMap, double scale)
   : Coefficient()
{
   // make a copy of the magic attribute-value map for later use
   materialMap = new std::map<int, double>(inputMap);
   scaleFactor = scale;
}

MeshDependentCoefficient::MeshDependentCoefficient(
   const MeshDependentCoefficient &cloneMe)
   : Coefficient()
{
   // make a copy of the magic attribute-value map for later use
   materialMap = new std::map<int, double>(*(cloneMe.materialMap));
   scaleFactor = cloneMe.scaleFactor;
}

double MeshDependentCoefficient::Eval(ElementTransformation &T,
                                      const IntegrationPoint &ip)
{
   // given the attribute, extract the coefficient value from the map
   std::map<int, double>::iterator it;
   int thisAtt = T.Attribute;
   double value;
   it = materialMap->find(thisAtt);
   if (it != materialMap->end())
   {
      value = it->second;
   }
   else
   {
      value = 0.0; // avoid compile warning
      std::cerr << "MeshDependentCoefficient attribute " << thisAtt
                << " not found" << std::endl;
      mfem_error();
   }

   return value*scaleFactor;
}

ScaledGFCoefficient::ScaledGFCoefficient(GridFunction *gf,
                                         MeshDependentCoefficient &input_mdc)
   : GridFunctionCoefficient(gf), mdc(input_mdc) {}

double ScaledGFCoefficient::Eval(ElementTransformation &T,
                                 const IntegrationPoint &ip)
{
   return mdc.Eval(T,ip) * GridFunctionCoefficient::Eval(T,ip);
}
double V_BCFunc_in(const Vector &x)
{
double phi = 10;
return phi;
}

double V_BCFunc_out(const Vector &x)
{
double phi = 0;
return phi;  
}
} // namespace electromagnetics

} // namespace mfem

#endif // MFEM_USE_MPI
