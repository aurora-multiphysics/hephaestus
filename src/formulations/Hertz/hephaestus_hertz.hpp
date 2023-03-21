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
#pragma once
#include <fstream>
#include <iostream>

#include "hertz_solver.hpp"
#include "inputs.hpp"
#include "postprocessors.hpp"

using namespace std;
using namespace mfem;
using namespace mfem::electromagnetics;

// Permittivity Functions
Coefficient *SetupPermittivityCoefficient();

static Vector pw_eps_(0);    // Piecewise permittivity values
static Vector ds_params_(0); // Center, Radius, and Permittivity
//                               of dielectric sphere
double dielectric_sphere(const Vector &);

// Permeability Function
Coefficient *SetupInvPermeabilityCoefficient();

static Vector pw_mu_(0);     // Piecewise permeability values
static Vector pw_mu_inv_(0); // Piecewise inverse permeability values
static Vector ms_params_(0); // Center, Inner and Outer Radii, and
//                               Permeability of magnetic shell
double magnetic_shell(const Vector &);
double magnetic_shell_inv(const Vector &x) { return 1.0 / magnetic_shell(x); }

// Conductivity Functions
Coefficient *SetupConductivityCoefficient();

static Vector pw_sigma_(0);  // Piecewise conductivity values
static Vector cs_params_(0); // Center, Radius, and Conductivity
//                               of conductive sphere
double conductive_sphere(const Vector &);

// Impedance
Coefficient *SetupAdmittanceCoefficient(const Mesh &mesh,
                                        const Array<int> &abcs);

static Vector pw_eta_(0);     // Piecewise impedance values
static Vector pw_eta_inv_(0); // Piecewise inverse impedance values

// Current Density Function
static Vector do_params_(0); // Axis Start, Axis End, Rod Radius,
//                               Total Current of Rod
void dipole_oscillator(const Vector &x, Vector &j);
void j_src(const Vector &x, Vector &j) { dipole_oscillator(x, j); }

// Electric Field Boundary Condition: The following function returns zero but
// any function could be used.
// void e_bc_r(const Vector &x, Vector &E);
// void e_bc_i(const Vector &x, Vector &E);

static double freq_ = 9.3e9; // 10/2pi

// double port_length_vector[3] = {24.76e-2, 0.0, 0.0};
// double port_width_vector[3] = {0.0, 12.38e-2, 0.0};
double port_length_vector[3] = {0.0, 22.86e-3, 0.0};
double port_width_vector[3] = {0.0, 0.0, 10.16e-3};

Vector a1Vec(port_length_vector, 3);
Vector a2Vec(port_width_vector, 3);
Vector a3Vec;
Vector a2xa3;
Vector a3xa1;
double kc;
double k0;
complex<double> k_;

void RWTE10(const Vector &x, vector<complex<double>> &E) {
  complex<double> zi = complex<double>(0., 1.);
  double port_length_vector[3] = {0.0, 22.86e-3, 0.0};
  double port_width_vector[3] = {0.0, 0.0, 10.16e-3};
  double omega_ = 2 * 3.14159 * freq_;
  Vector a1Vec(port_length_vector, 3);
  Vector a2Vec(port_width_vector, 3);
  Vector a3Vec;
  Vector a2xa3;
  Vector a3xa1;

  hephaestus::cross_product(a1Vec, a2Vec, a3Vec);
  hephaestus::cross_product(a2Vec, a3Vec, a2xa3);
  hephaestus::cross_product(a3Vec, a1Vec, a3xa1);

  Vector k_a = a2xa3;
  Vector k_c = a3Vec;
  double V = InnerProduct(a1Vec, a2xa3);
  k_a *= M_PI / V;
  k_c *= k_.imag() / a3Vec.Norml2();

  Vector E_hat;
  hephaestus::cross_product(k_c, k_a, E_hat);
  E_hat *= 1.0 / E_hat.Norml2();

  double E0(
      sqrt(2 * omega_ * mu0_ / (a1Vec.Norml2() * a2Vec.Norml2() * k_.imag())));
  complex<double> E_mag =
      E0 * sin(InnerProduct(k_a, x)) * exp(-zi * InnerProduct(k_c, x));

  E[0] = E_mag * E_hat(1);
  E[1] = E_mag * E_hat(2);
  E[2] = E_mag * E_hat(0);
}

void RWTE10_real(const Vector &x, Vector &v) {
  vector<complex<double>> Eval(x.Size());
  RWTE10(x, Eval);
  for (int i = 0; i < x.Size(); ++i) {
    v(i) = -2 * k_.imag() * Eval[i].imag() / mu0_;
  }
}
void RWTE10_imag(const Vector &x, Vector &v) {
  vector<complex<double>> Eval(x.Size());
  RWTE10(x, Eval);
  for (int i = 0; i < x.Size(); ++i) {
    v(i) = 2 * k_.imag() * Eval[i].real() / mu0_;
  }
}

int hertz_solve(int argc, char *argv[], hephaestus::Inputs inputs);

// void e_bc_r(const Vector &x, Vector &E)
// {
//    E.SetSize(3);
//    E = 0.0;
// }

// void e_bc_i(const Vector &x, Vector &E)
// {
//    E.SetSize(3);
//    E = 0.0;
// }
