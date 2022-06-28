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
//  mpirun -np 4 valgrind Joint -m  Joint_slab2.msh
// lookup-table-based nonlinear modelling produced, temperature-dependent electrical conductivity
#include "Joint_solver.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;
using namespace mfem::electromagnetics;

double muInv(const Vector & x) {const double mu0_ = 1.2566e-6; return 1.0/mu0_; }
double muFunc(const Vector & x) {const double mu0_ = 1.2566e-6; return mu0_; }
double T_sigma(double & x) {const double sigma = 1.0e-10; return sigma; }
static double mj_ = 0.0;
static double sj_ = 0.0;
static double wj_ = 0.0;
double ramp_current1(double t)
{     
   return  500 * t;// current is ramped up with time from zero 10A/s
}
double ItrFunc( double t) { 
   double I = (t < 2) ? ramp_current1(t) : 750;
   return I ; }

double T_sigma_copper( double T) { 
   const double RRR = 100; 
   double rho0 = 15.5e-9/RRR;
   double rho_i = (1.171e-17*pow(T, 4.49))/ (1+1.171e-17*3.841e10*pow(T, 3.35)*exp(-pow( 50/T, 6.428)));
   double rho_io = 0.4531*rho_i*rho0/(rho0 + rho_i) ;
   double sigma = 1/(rho0 + rho_i + rho_io);
   return sigma ; }
double T_sigma_HTS( double T) { 
   double sigma = 1e16;
   return sigma ; } 

double T_sigma_X( double T) { 
   double sigma = 1.0 ;
   return sigma  ; }      

void dBdt_xyz_BC(const Vector &x, double t, Vector&dBdt)
{
    dBdt.SetSize(3);
    if (x[0] == 0 || x[0] == 2)
    {
     dBdt(0) = 0.0*0.5/0.5;// dt = 0.5
     dBdt(1) = 0.0;
     dBdt(2) = 0.0;
      //dBdt =0.0;
    }
    else if (x[1] == 0 || x[1] == 2)
    {
     dBdt =0.0;
    }
    else if (x[2] == 0  || x[2] == 2)
    {
     dBdt =0.0;  
    }
    else
    {dBdt = 0.0;}
   
}
int electromagnetics::SOLVER_PRINT_LEVEL = 0;
int electromagnetics::STATIC_COND        = 0;

int main(int argc, char *argv[])
{
// 1. Initialize MPI.
   int num_procs, myid;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   
// Parse command-line options.
   const char *mesh_file = "../torus2.msh";
   int Order = 1;
   int serial_ref_levels = 1;
   int parallel_ref_levels = 0;
   bool visit = true;

   double dtsf = 0.95;
   double ti = 0.0;
   double ts = 0.5;
   double tf =  4;
   const char *petscrc_file = "";
   int ode_solver_type = 1;
   double t_final = 100.0;
   double dt = 0.5;
   double mu0 = 1.257e-6;
   double sigma = 1.0e10;
   double Tcapacity = 1.0;
   double Tconductivity = 0.01;
   double freq = 1.0/60.0;
   const char *basename = "Joint EM-thermal";
   OptionsParser args(argc, argv);

   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
  // args.AddOption(&petscrc_file, "-petscopts", "--petscopts",
               //   "PetscOptions file to use.");

   args.Parse();
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(cout);
   }
   // initilize PETSC
  // MFEMInitializePetsc(NULL,NULL,petscrc_file,NULL);

   // Read the (serial) mesh from the given mesh file on all processors.  We can
   // handle triangular, quadrilateral, tetrahedral, hexahedral, surface and
   // volume meshes with the same code. 
   double mu = mu0  ;
   mj_  = mu;
   sj_  = sigma;
   wj_  = 2.0*M_PI*freq;

   if (myid == 1)
   {
      cout << "\nSkin depth sqrt(2.0/(wj*mj*sj)) = " << sqrt(2.0/(wj_*mj_*sj_))
           << "\nSkin depth sqrt(2.0*dt/(mj*sj)) = " << sqrt(2.0*dt/(mj_*sj_))
           << endl;
   }

   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();
   if (myid == 0){ cout << " mesh used " << dim << endl; }
   Array<int> ess_bdr(mesh->bdr_attributes.Max());
   Array<int> domains_(mesh->attributes.Max());
   Array<int> thermal_ess_bdr(mesh->bdr_attributes.Max());
   Array<int> poisson_ess_bdr(mesh->bdr_attributes.Max());
   Array<int> poisson_ess_bdr_in(mesh->bdr_attributes.Max()),poisson_ess_bdr_out(mesh->bdr_attributes.Max()),poisson_ess_bdr_outer(mesh->bdr_attributes.Max());
   Array<int> magnetic_nat_bdr_x(mesh->bdr_attributes.Max()),magnetic_nat_bdr_y(mesh->bdr_attributes.Max()), magnetic_nat_bdr_z(mesh->bdr_attributes.Max()) ;
  
      ess_bdr = 1;

      // Same as above, but this is for the thermal operator for HDiv
      // formulation the essential BC is the flux
      poisson_ess_bdr = 1;
      poisson_ess_bdr[2] = 0; // boundary attribute 1 (index 0) is fixed
      
      poisson_ess_bdr_in = 0;

      poisson_ess_bdr_in[4] = 1;

      poisson_ess_bdr_out = 0;
      poisson_ess_bdr_out[5] = 1;

      poisson_ess_bdr_outer = 0;
      poisson_ess_bdr_outer[0] = 1;
      poisson_ess_bdr_outer[1] = 1;
      poisson_ess_bdr_outer[3] = 1;
      
      magnetic_nat_bdr_x = 0; // perpendicular magnetic flux B to the outer surfaces, in x direction 
      magnetic_nat_bdr_x[0] = 1;

      magnetic_nat_bdr_y = 1;
      magnetic_nat_bdr_y[0] = 0;
      magnetic_nat_bdr_y[3] = 0;
      
      magnetic_nat_bdr_z = 0;
      magnetic_nat_bdr_z[3] = 1;

      thermal_ess_bdr = 1; // 

      
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   pmesh->ReorientTetMesh();

   // 11. Rebalance the mesh. Since the mesh was adaptively refined in a
   //     non-uniform way it will be computationally unbalanced.
   if (pmesh->Nonconforming())
   {
      pmesh->Rebalance();
   }
   if (myid == 0){ cout << " number of boundaries  = " << pmesh->bdr_attributes.Max() << endl; }
   if (myid == 0){ cout << " number of domains  = " << pmesh->attributes.Max() << endl; }
   
    /* the ampersand is actually optional */

   std::map<int, double> sigmaMap, InvTcondMap, TcapMap, InvTcapMap;
   double sigma_Air, sigma_copper, sigma_aluminium,sigma_X ;
   double Tcond_Air,Tcond_Copper,Tcond_Aluminium, Tcond_X ;
   double Tcap_Air,Tcap_Copper, Tcap_Aluminium,Tcap_X ;
   
      sigma_Air     = 1.0e-10 * sigma;
      Tcond_Air     = 1.0e-3  * Tconductivity;
      Tcap_Air      = 10000    * Tcapacity;
      sigma_copper = 0.1*sigma;
      sigma_aluminium = 0.01*sigma;
      sigma_X = 100 *sigma;

      Tcond_Copper = Tconductivity;
      Tcond_Aluminium = 0.1*Tconductivity;
      Tcond_X = 0.01*Tconductivity;
      
      Tcap_Copper = Tcapacity*0.385;
      Tcap_Aluminium = 0.89*Tcapacity;
      Tcap_X = 2 *Tcapacity;
     
      sigmaMap.insert(pair<int, double>(1, sigma_Air));
      sigmaMap.insert(pair<int, double>(2, sigma_X));
      sigmaMap.insert(pair<int, double>(3, sigma_X));
      sigmaMap.insert(pair<int, double>(4, sigma_copper));
      sigmaMap.insert(pair<int, double>(5, sigma_copper));
      //SigmaT.insert(pair<int, double*>(1, *SigmaTT);
      //sigmaMap.insert(pair<int, double>(5, sigma_aluminium));

      InvTcondMap.insert(pair<int, double>(1, 1.0/Tcond_Air));
      InvTcondMap.insert(pair<int, double>(2, 1.0/Tcond_X));
      InvTcondMap.insert(pair<int, double>(3, 1.0/Tcond_X));
      InvTcondMap.insert(pair<int, double>(4, 1.0/Tcond_Copper));
      InvTcondMap.insert(pair<int, double>(5, 1.0/Tcond_Copper));
      //InvTcondMap.insert(pair<int, double>(5, 1.0/Tcond_Aluminium));

      TcapMap.insert(pair<int, double>(1, Tcap_Air));
      TcapMap.insert(pair<int, double>(2, Tcap_X));
      TcapMap.insert(pair<int, double>(3, Tcap_X));
      TcapMap.insert(pair<int, double>(4, Tcap_Copper));
      TcapMap.insert(pair<int, double>(5, Tcap_Copper));

      InvTcapMap.insert(pair<int, double>(1, 1.0/Tcap_Air));
      InvTcapMap.insert(pair<int, double>(2, 1.0/Tcap_X));
      InvTcapMap.insert(pair<int, double>(3, 1.0/Tcap_X));
      InvTcapMap.insert(pair<int, double>(4, 1.0/Tcap_Copper));
      InvTcapMap.insert(pair<int, double>(5, 1.0/Tcap_Copper));
   FiniteElementCollection *fec_ND; 
   FiniteElementCollection *fec_RT; 
   FiniteElementCollection *fec_L2;
   FiniteElementCollection *fec_H1;
   fec_ND= new ND_FECollection(Order, pmesh->Dimension());
   fec_RT= new RT_FECollection(Order, pmesh->Dimension());
   fec_L2= new L2_FECollection(Order, pmesh->Dimension());
   fec_H1= new H1_FECollection(Order, pmesh->Dimension());
   ParFiniteElementSpace HCurlFESpace(pmesh,fec_ND);
   ParFiniteElementSpace HDivFESpace(pmesh,fec_RT);
   ParFiniteElementSpace L2FESpace(pmesh, fec_L2);
   ParFiniteElementSpace HGradFESpace(pmesh, fec_H1);

   HYPRE_BigInt glob_size_l2 = L2FESpace.GlobalTrueVSize();
   HYPRE_BigInt glob_size_nd = HCurlFESpace.GlobalTrueVSize();
   HYPRE_BigInt glob_size_rt = HDivFESpace.GlobalTrueVSize();
   HYPRE_BigInt glob_size_h1 = HGradFESpace.GlobalTrueVSize();

   if (myid == 0)
   {
      cout << "Number of Temperature Flux unknowns:  " << glob_size_rt << endl;
      cout << "Number of Temperature unknowns:       " << glob_size_l2 << endl;
      cout << "Number of Electric Field unknowns:    " << glob_size_nd << endl;
      cout << "Number of Magnetic Field unknowns:    " << glob_size_rt << endl;
      cout << "Number of Electrostatic unknowns:     " << glob_size_h1 << endl;
   }

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
   // The BlockVector is a large contiguous chunk of memory for storing required
   // data for the hypre vectors, in this case: the temperature L2, the T-flux
   // HDiv, the E-field HCurl, and the B-field HDiv, and scalar potential P.
   BlockVector F(true_offset);

   // grid functions E, B, T, F, P, and w which is the Joule heating
   ParGridFunction P_gf, E_gf, B_gf, T_gf, F_gf, w_gf, J_gf, J_in;
   P_gf.MakeRef(&HGradFESpace,F,   true_offset[0]);
   T_gf.MakeRef(&L2FESpace,F,   true_offset[1]);
   F_gf.MakeRef(&HDivFESpace,F, true_offset[2]);
   J_gf.MakeRef(&HDivFESpace,F, true_offset[3]);
   E_gf.MakeRef(&HCurlFESpace,F,true_offset[4]);
   B_gf.MakeRef(&HDivFESpace,F, true_offset[5]);
   w_gf.MakeRef(&L2FESpace,F,   true_offset[6]);
   J_in.MakeRef(&HDivFESpace,F, true_offset[7]);

   int size_Hcurl =  HCurlFESpace.GetTrueVSize();
    if ( myid == 0 )  {  cout << "Number of local H(Curl) unknowns: " << size_Hcurl << endl;  }
   HYPRE_Int size_rt = HDivFESpace.GlobalTrueVSize();
    if ( myid == 0 )  {  cout << "Number of global H(Div) unknowns: " << size_rt << endl;   }
    if ( myid == 0 )  {  cout << "Number of local L2 unknowns: " << L2FESpace.GetTrueVSize() << endl;   }
   // Create the Electromagnetic solver
   
   Joint_solver oper(true_offset[8], L2FESpace, HCurlFESpace, HDivFESpace, HGradFESpace, domains_, 
                                   ess_bdr, thermal_ess_bdr, poisson_ess_bdr,poisson_ess_bdr_in, poisson_ess_bdr_out, poisson_ess_bdr_outer,
                                   magnetic_nat_bdr_x, magnetic_nat_bdr_y, magnetic_nat_bdr_z,
                                   mu, sigmaMap, TcapMap, InvTcapMap, InvTcondMap, ItrFunc, dBdt_xyz_BC, T_sigma_copper,T_sigma_HTS, T_sigma_X );
    
    if ( myid == 0 )  {  cout << "Starting initialization Magnet solver." << endl;  }
   oper.Init(F);
   // Set the largest stable time step
   double dtmax = 0.001 ;

   // Create the ODE solver
   BackwardEulerSolver BESolver;
   BESolver.Init(oper);
    if ( myid == 0 )  {  cout << "Initialization ODE solver finished." << endl; }
   // Initialize VisIt visualization
   double t = ti;
   oper.SetTime(t);

   
   // Write initial fields to disk for VisIt
   VisItDataCollection visit_dc(basename, pmesh);
      
      visit_dc.RegisterField("B", &B_gf);
      visit_dc.RegisterField("T", &T_gf);
      visit_dc.RegisterField("w", &w_gf);
      visit_dc.RegisterField("J", &J_gf);
      visit_dc.RegisterField("J_induced", &J_in);
      visit_dc.SetCycle(0);
      visit_dc.SetTime(0.0);
      visit_dc.Save();
   
   // The main time evolution loop.  
   int it = 0;
   t = 0;
   while (t < tf)
   {
      // Run the simulation until a snapshot is needed
       BESolver.Step(F, t, dt);  // Step() includes t += dt     H = H +dHdt*dt

      // Update local DoFs with current true DoFs
      // MagnetEM.SyncGridFuncs();
      visit_dc.SetCycle(it);
      visit_dc.SetTime(t);
      visit_dc.Save();
      // Write fields to disk for VisIt
   it++;  
   }

   delete fec_H1;
   delete fec_L2;
   delete fec_ND;
   delete fec_RT;
   delete pmesh;
 // MFEMFinalizePetsc();
  //MPI_Finalize();
   return 0;
}