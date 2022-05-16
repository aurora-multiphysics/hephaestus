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
//
//   -----------------------------------------------------------------------
//   Hertz Miniapp:  Simple Frequency-Domain Electromagnetic Simulation Code
//   -----------------------------------------------------------------------
//
//   Assumes that all sources and boundary conditions oscillate with the same
//   frequency although not necessarily in phase with one another.  This
//   assumptions implies that we can factor out the time dependence which we
//   take to be of the form exp(i omega t).  With these assumptions we can
//   write the Maxwell equations in the form:
//
//   i omega epsilon E = Curl mu^{-1} B - J - sigma E
//   i omega B         = - Curl E
//
//   Which combine to yield:
//
//   Curl mu^{-1} Curl E - omega^2 epsilon E + i omega sigma E = - i omega J
//
//   We discretize this equation with H(Curl) a.k.a Nedelec basis
//   functions.  The curl curl operator must be handled with
//   integration by parts which yields a surface integral:
//
//   (W, Curl mu^{-1} Curl E) = (Curl W, mu^{-1} Curl E)
//               + (W, n x (mu^{-1} Curl E))_{\Gamma}
//
//   or
//
//   (W, Curl mu^{-1} Curl E) = (Curl W, mu^{-1} Curl E)
//               - i omega (W, n x H)_{\Gamma}
//
//   For plane waves
//     omega B = - k x E
//     omega D = k x H, assuming n x k = 0 => n x H = omega epsilon E / |k|
//
//   c = omega/|k|
//
//   (W, Curl mu^{-1} Curl E) = (Curl W, mu^{-1} Curl E)
//               - i omega sqrt{epsilon/mu} (W, E)_{\Gamma}
//
//
// Compile with: make hertz
//
// Sample runs:
//
//   By default the sources and fields are all zero
//     mpirun -np 4 hertz
//
//   Current source in a metal sphere
//     mpirun -np 4 hertz -m ../../data/ball-nurbs.mesh -rs 2
//                        -dbcs '-1' -f 3e8 -herm
//                        -do '-0.3 0.0 0.0 0.3 0.0 0.0 0.1 1 .5 .5'
//
//   Current source in a sphere with absorbing boundary conditions
//     mpirun -np 4 hertz -m ../../data/ball-nurbs.mesh -rs 2
//                        -abcs '-1' -f 3e8
//                        -do '-0.3 0.0 0.0 0.3 0.0 0.0 0.1 1 .5 .5'
//
//   Current source in a metal sphere with dielectric and conducting materials
//     mpirun -np 4 hertz -m ../../data/ball-nurbs.mesh -rs 2
//                        -dbcs '-1' -f 3e8
//                        -do '-0.3 0.0 0.0 0.3 0.0 0.0 0.1 1 .5 .5'
//                        -cs '0.0 0.0 -0.5 .2 10'
//                        -ds '0.0 0.0 0.5 .2 10'
//
//   Current source in a metal box
//     mpirun -np 4 hertz -m ../../data/fichera.mesh -rs 3
//                        -dbcs '-1' -f 3e8
//                        -do '-0.5 -0.5 0.0 -0.5 -0.5 1.0 0.1 1 .5 1'
//
//   Current source with a mixture of absorbing and reflecting boundaries
//     mpirun -np 4 hertz -m ../../data/fichera.mesh -rs 3
//                        -do '-0.5 -0.5 0.0 -0.5 -0.5 1.0 0.1 1 .5 1'
//                        -dbcs '4 8 19 21' -abcs '5 18' -f 3e8
//
#include "hephaestus_hertz.hpp"

using namespace std;
using namespace mfem;
using namespace mfem::electromagnetics;

int hertz_solve(int argc, char *argv[], hephaestus::Inputs inputs) {
  int world_rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Parse command-line options.
  const char *formulation = "Hertz";
  const char *mesh_file = inputs.mesh_file.c_str();
  int order = 1;
  int maxit = 1;
  int serial_ref_levels = 0;
  int parallel_ref_levels = 0;
  int sol = 6;
  int prec = 4;
  bool herm_conv = true;
  bool visualization = true;
  bool visit = true;

  Array<int> abcs;
  Array<int> dbcs;

  SolverOptions solOpts;
  solOpts.maxIter = 1000;
  solOpts.kDim = 50;
  solOpts.printLvl = 1;
  solOpts.relTol = 1e-2;
  solOpts.euLvl = 2;

  OptionsParser args(argc, argv);
  args.AddOption(&formulation, "-form", "--formulation",
                 "Name of formulation to use during solve.");
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&serial_ref_levels, "-rs", "--serial-ref-levels",
                 "Number of serial refinement levels.");
  args.AddOption(&parallel_ref_levels, "-rp", "--parallel-ref-levels",
                 "Number of parallel refinement levels.");
  args.AddOption(&freq_, "-f", "--frequency",
                 "Frequency in Hertz (of course...)");
  args.AddOption(&prec, "-pc", "--precond",
                 "Preconditioner: 1 - Diagonal Scaling, 2 - ParaSails, "
                 "3 - Euclid, 4 - AMS");
  args.AddOption(&sol, "-s", "--solver",
                 "Solver: 1 - GMRES, 2 - FGMRES, 3 - MINRES"
#ifdef MFEM_USE_SUPERLU
                 ", 4 - SuperLU"
#endif
#ifdef MFEM_USE_STRUMPACK
                 ", 5 - STRUMPACK"
#endif
#ifdef MFEM_USE_MUMPS
                 ", 6 - MUMPS"
#endif
  );
  args.AddOption(&solOpts.maxIter, "-sol-it", "--solver-iterations",
                 "Maximum number of solver iterations.");
  args.AddOption(&solOpts.kDim, "-sol-k-dim", "--solver-krylov-dimension",
                 "Krylov space dimension for GMRES and FGMRES.");
  args.AddOption(&solOpts.relTol, "-sol-tol", "--solver-tolerance",
                 "Relative tolerance for GMRES or FGMRES.");
  args.AddOption(&solOpts.printLvl, "-sol-prnt-lvl", "--solver-print-level",
                 "Logging level for solvers.");
  args.AddOption(&solOpts.euLvl, "-eu-lvl", "--euclid-level",
                 "Euclid factorization level for ILU(k).");
  args.AddOption(&pw_eps_, "-pwe", "--piecewise-eps",
                 "Piecewise values of Permittivity");
  args.AddOption(&ds_params_, "-ds", "--dielectric-sphere-params",
                 "Center, Radius, and Permittivity of Dielectric Sphere");
  args.AddOption(&pw_mu_, "-pwm", "--piecewise-mu",
                 "Piecewise values of Permeability");
  args.AddOption(&ms_params_, "-ms", "--magnetic-shell-params",
                 "Center, Inner Radius, Outer Radius, "
                 "and Permeability of Magnetic Shell");
  args.AddOption(&pw_sigma_, "-pws", "--piecewise-sigma",
                 "Piecewise values of Conductivity");
  args.AddOption(&cs_params_, "-cs", "--conductive-sphere-params",
                 "Center, Radius, and Conductivity of Conductive Sphere");
  args.AddOption(&pw_eta_, "-pwz", "--piecewise-eta",
                 "Piecewise values of Impedance (one value per abc surface)");
  args.AddOption(&do_params_, "-do", "--dipole-oscillator-params",
                 "Axis End Points, Radius, and Amplitude");
  args.AddOption(&abcs, "-abcs", "--absorbing-bc-surf",
                 "Absorbing Boundary Condition Surfaces");
  args.AddOption(&dbcs, "-dbcs", "--dirichlet-bc-surf",
                 "Dirichlet Boundary Condition Surfaces");
  args.AddOption(&maxit, "-maxit", "--max-amr-iterations",
                 "Max number of iterations in the main AMR loop.");
  args.AddOption(&herm_conv, "-herm", "--hermitian", "-no-herm",
                 "--no-hermitian", "Use convention for Hermitian operators.");
  args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                 "--no-visualization",
                 "Enable or disable GLVis visualization.");
  args.AddOption(&visit, "-visit", "--visit", "-no-visit", "--no-visit",
                 "Enable or disable VisIt visualization.");
  args.Parse();
  if (!args.Good()) {
    if (world_rank == 0) {
      args.PrintUsage(cout);
    }
    return 1;
  }
  if (world_rank == 0) {
    args.PrintOptions(cout);
  }

  ComplexOperator::Convention conv =
      herm_conv ? ComplexOperator::HERMITIAN : ComplexOperator::BLOCK_SYMMETRIC;

  // Read the (serial) mesh from the given mesh file on all processors.  We
  // can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
  // and volume meshes with the same code.
  Mesh *mesh = new Mesh(mesh_file, 1, 1);
  if (world_rank == 0) {
    cout << "Starting initialization." << endl;
  }

  // Project a NURBS mesh to a piecewise-quadratic curved mesh
  if (mesh->NURBSext) {
    mesh->UniformRefinement();
    if (serial_ref_levels > 0) {
      serial_ref_levels--;
    }

    mesh->SetCurvature(2);
  }

  // Ensure that quad and hex meshes are treated as non-conforming.
  mesh->EnsureNCMesh();

  // Refine the serial mesh on all processors to increase the resolution. In
  // this example we do 'ref_levels' of uniform refinement.
  for (int l = 0; l < serial_ref_levels; l++) {
    mesh->UniformRefinement();
  }

  // Define a parallel mesh by a partitioning of the serial mesh. Refine
  // this mesh further in parallel to increase the resolution. Once the
  // parallel mesh is defined, the serial mesh can be deleted.
  ParMesh pmesh(MPI_COMM_WORLD, *mesh);
  delete mesh;

  // Refine this mesh in parallel to increase the resolution.
  int par_ref_levels = parallel_ref_levels;
  for (int l = 0; l < par_ref_levels; l++) {
    pmesh.UniformRefinement();
  }

  // Create a coefficient describing the dielectric permittivity
  Coefficient *epsCoef = SetupPermittivityCoefficient();

  // Create a coefficient describing the magnetic permeability
  Coefficient *muInvCoef = SetupInvPermeabilityCoefficient();

  // Create a coefficient describing the electrical conductivity
  Coefficient *sigmaCoef = SetupConductivityCoefficient();

  // Create a coefficient describing the surface admittance
  // Coefficient * etaInvCoef = SetupAdmittanceCoefficient(pmesh, abcs);

  // Coefficient * etaInvCoef = new PWConstCoefficient(pw_eta_inv_);;

  // General Robin: a*(nxcurlE) + b*nxnxE = U
  // Free space a = 1, b = -i omega sqrt(eps0*mu0), U =
  kc = M_PI / a1Vec.Norml2();
  freq_ = 9.3e9;
  double omega_ = 2.0 * M_PI * freq_;
  k0 = omega_ * sqrt(epsilon0_ * mu0_);
  k_ = complex<double>(0., sqrt(k0 * k0 - kc * kc));

  // Robin coefficient, but already handled here.
  Coefficient *etaInvCoef =
      new ConstantCoefficient(-k_.imag() / (mu0_ * omega_));

  hephaestus::BCMap bc_map = inputs.bc_map;

  hephaestus::BoundaryCondition waveguide_ports(std::string("robin_1"),
                                                Array<int>({1, 2}));
  Array<int> wgi_in_attr(1);
  wgi_in_attr[0] = 1;
  hephaestus::IntegratedBC waveguide_in(std::string("wgi"), wgi_in_attr);
  VectorFunctionCoefficient UReal(pmesh.SpaceDimension(), RWTE10_real);
  VectorFunctionCoefficient UImag(pmesh.SpaceDimension(), RWTE10_imag);
  waveguide_in.lfi_re = new mfem::VectorFEBoundaryTangentLFIntegrator(UReal);
  waveguide_in.lfi_im = new mfem::VectorFEBoundaryTangentLFIntegrator(UImag);
  waveguide_in.markers.SetSize(pmesh.bdr_attributes.Max());
  waveguide_in.markers = 0;
  waveguide_in.markers[1] = 1;

  dbcs.SetSize(1);
  dbcs = 0;
  dbcs[0] = 1;

  abcs = Array<int>({2, 3});

  bc_map["Neumann"] = new hephaestus::IntegratedBC(waveguide_in);

  hephaestus::VectorFunctionDirichletBC *tangential_E_bc =
      dynamic_cast<hephaestus::VectorFunctionDirichletBC *>(
          bc_map["tangential_E"]);

  // Create the Magnetostatic solver
  HertzSolver Hertz(pmesh, order, freq_, (HertzSolver::SolverType)sol, solOpts,
                    (HertzSolver::PrecondType)prec, conv, *epsCoef, *muInvCoef,
                    sigmaCoef, etaInvCoef, bc_map, abcs, dbcs,
                    (do_params_.Size() > 0) ? j_src : NULL, NULL);

  // Initialize GLVis visualization
  if (visualization) {
    Hertz.InitializeGLVis();
  }

  // Initialize VisIt visualization
  VisItDataCollection visit_dc("Hertz-AMR-Parallel", &pmesh);

  if (visit) {
    Hertz.RegisterVisItFields(visit_dc);
  }
  if (world_rank == 0) {
    cout << "Initialization done." << endl;
  }

  // The main AMR loop. In each iteration we solve the problem on the current
  // mesh, visualize the solution, estimate the error on all elements, refine
  // the worst elements and update all objects to work with the new mesh. We
  // refine until the maximum number of dofs in the Nedelec finite element
  // space reaches 10 million.
  const int max_dofs = 10000000;
  for (int it = 1; it <= maxit; it++) {
    if (world_rank == 0) {
      cout << "\nAMR Iteration " << it << endl;
    }

    // Display the current number of DoFs in each finite element space
    Hertz.PrintSizes();

    // Assemble all forms
    Hertz.Assemble();

    // Solve the system and compute any auxiliary fields
    Hertz.Solve();

    // Determine the current size of the linear system
    int prob_size = Hertz.GetProblemSize();

    // Write fields to disk for VisIt
    if (visit) {
      Hertz.WriteVisItFields(it);
    }

    // Send the solution by socket to a GLVis server.
    if (visualization) {
      Hertz.DisplayToGLVis();
    }

    if (world_rank == 0) {
      cout << "AMR iteration " << it << " complete." << endl;
    }

    // Check stopping criteria
    if (prob_size > max_dofs) {
      if (world_rank == 0) {
        cout << "Reached maximum number of dofs, exiting..." << endl;
      }
      break;
    }
    if (it == maxit) {
      break;
    }

    // Wait for user input. Ask every 10th iteration.
    char c = 'c';
    if (world_rank == 0 && (it % 10 == 0)) {
      cout << "press (q)uit or (c)ontinue --> " << flush;
      cin >> c;
    }
    MPI_Bcast(&c, 1, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (c != 'c') {
      break;
    }

    // Estimate element errors using the Zienkiewicz-Zhu error estimator.
    Vector errors(pmesh.GetNE());
    Hertz.GetErrorEstimates(errors);

    double local_max_err = errors.Max();
    double global_max_err;
    MPI_Allreduce(&local_max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX,
                  pmesh.GetComm());

    // Refine the elements whose error is larger than a fraction of the
    // maximum element error.
    const double frac = 0.5;
    double threshold = frac * global_max_err;
    if (world_rank == 0) {
      cout << "Refining ..." << endl;
    }
    pmesh.RefineByError(errors, threshold);

    // Update the magnetostatic solver to reflect the new state of the mesh.
    Hertz.Update();

    if (pmesh.Nonconforming() && world_size > 1 && false) {
      if (world_rank == 0) {
        cout << "Rebalancing ..." << endl;
      }
      pmesh.Rebalance();

      // Update again after rebalancing
      Hertz.Update();
    }
  }

  // Send the solution by socket to a GLVis server.
  if (visualization) {
    Hertz.DisplayAnimationToGLVis();
  }

  delete epsCoef;
  delete muInvCoef;
  delete sigmaCoef;
  delete etaInvCoef;

  return 0;
}

// The Permittivity is a required coefficient which may be defined in
// various ways so we'll determine the appropriate coefficient type here.
Coefficient *SetupPermittivityCoefficient() {
  Coefficient *coef = NULL;

  if (ds_params_.Size() > 0) {
    coef = new FunctionCoefficient(dielectric_sphere);
  } else if (pw_eps_.Size() > 0) {
    coef = new PWConstCoefficient(pw_eps_);
  } else {
    coef = new ConstantCoefficient(epsilon0_);
  }

  return coef;
}

// The Permeability is a required coefficient which may be defined in
// various ways so we'll determine the appropriate coefficient type here.
Coefficient *SetupInvPermeabilityCoefficient() {
  Coefficient *coef = NULL;

  if (ms_params_.Size() > 0) {
    coef = new FunctionCoefficient(magnetic_shell_inv);
  } else if (pw_mu_.Size() > 0) {
    pw_mu_inv_.SetSize(pw_mu_.Size());
    for (int i = 0; i < pw_mu_.Size(); i++) {
      MFEM_ASSERT(pw_mu_[i] > 0.0, "permeability values must be positive");
      pw_mu_inv_[i] = 1.0 / pw_mu_[i];
    }
    coef = new PWConstCoefficient(pw_mu_inv_);
  } else {
    coef = new ConstantCoefficient(1.0 / mu0_);
  }

  return coef;
}

// The Conductivity is an optional coefficient which may be defined in
// various ways so we'll determine the appropriate coefficient type here.
Coefficient *SetupConductivityCoefficient() {
  Coefficient *coef = NULL;

  if (cs_params_.Size() > 0) {
    coef = new FunctionCoefficient(conductive_sphere);
  } else if (pw_sigma_.Size() > 0) {
    coef = new PWConstCoefficient(pw_sigma_);
  }

  return coef;
}

// The Admittance is an optional coefficient defined on boundary surfaces which
// can be used in conjunction with absorbing boundary conditions.
Coefficient *SetupAdmittanceCoefficient(const Mesh &mesh,
                                        const Array<int> &abcs) {
  Coefficient *coef = NULL;

  if (pw_eta_.Size() > 0) {
    MFEM_VERIFY(pw_eta_.Size() == abcs.Size(),
                "Each impedance value must be associated with exactly one "
                "absorbing boundary surface.");

    pw_eta_inv_.SetSize(mesh.bdr_attributes.Size());

    if (abcs[0] == -1) {
      pw_eta_inv_ = 1.0 / pw_eta_[0];
    } else {
      pw_eta_inv_ = 0.0;

      for (int i = 0; i < pw_eta_.Size(); i++) {
        pw_eta_inv_[abcs[i] - 1] = 1.0 / pw_eta_[i];
      }
    }
    coef = new PWConstCoefficient(pw_eta_inv_);
  }

  return coef;
}

// A sphere with constant permittivity.  The sphere has a radius,
// center, and permittivity specified on the command line and stored
// in ds_params_.
double dielectric_sphere(const Vector &x) {
  double r2 = 0.0;

  for (int i = 0; i < x.Size(); i++) {
    r2 += (x(i) - ds_params_(i)) * (x(i) - ds_params_(i));
  }

  if (sqrt(r2) <= ds_params_(x.Size())) {
    return ds_params_(x.Size() + 1) * epsilon0_;
  }
  return epsilon0_;
}

// A spherical shell with constant permeability.  The sphere has inner
// and outer radii, center, and relative permeability specified on the
// command line and stored in ms_params_.
double magnetic_shell(const Vector &x) {
  double r2 = 0.0;

  for (int i = 0; i < x.Size(); i++) {
    r2 += (x(i) - ms_params_(i)) * (x(i) - ms_params_(i));
  }

  if (sqrt(r2) >= ms_params_(x.Size()) &&
      sqrt(r2) <= ms_params_(x.Size() + 1)) {
    return mu0_ * ms_params_(x.Size() + 2);
  }
  return mu0_;
}

// A sphere with constant conductivity.  The sphere has a radius,
// center, and conductivity specified on the command line and stored
// in ls_params_.
double conductive_sphere(const Vector &x) {
  double r2 = 0.0;

  for (int i = 0; i < x.Size(); i++) {
    r2 += (x(i) - cs_params_(i)) * (x(i) - cs_params_(i));
  }

  if (sqrt(r2) <= cs_params_(x.Size())) {
    return cs_params_(x.Size() + 1);
  }
  return 0.0;
}

// A cylindrical rod of current density.  The rod has two axis end
// points, a radus, a current amplitude in Amperes.  All of these
// parameters are stored in do_params_.
void dipole_oscillator(const Vector &x, Vector &j) {
  MFEM_ASSERT(x.Size() == 3, "current source requires 3D space.");

  j.SetSize(x.Size());
  j = 0.0;

  Vector v(x.Size());  // Normalized Axis vector
  Vector xu(x.Size()); // x vector relative to the axis end-point

  xu = x;

  for (int i = 0; i < x.Size(); i++) {
    xu[i] -= do_params_[i];
    v[i] = do_params_[x.Size() + i] - do_params_[i];
  }

  double h = v.Norml2();

  if (h == 0.0) {
    return;
  }
  v /= h;

  double r = do_params_[2 * x.Size() + 0];
  double a = do_params_[2 * x.Size() + 1];

  double xv = xu * v;

  // Compute perpendicular vector from axis to x
  xu.Add(-xv, v);

  double xp = xu.Norml2();

  if (xv >= 0.0 && xv <= h && xp <= r) {
    j.Add(a, v);
  }
}
