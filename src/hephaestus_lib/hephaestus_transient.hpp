//            ---------------------------------------------------
//            ESolver:  Low-Frequency Electrodynamics Simulation
//            ---------------------------------------------------
//
// This miniapp solves low frequency magnetodynamics problems using the E
// formulation.
//
// (ν∇×E, ∇×E') - (σE, E') - (J0, E') - <(ν∇×E) × n, E'> = 0
// -(J0, ∇ V') + <n.J, V'> = 0

#pragma once
#include "../common/pfem_extras.hpp"
#include "eb_dual_solver.hpp"
#include "hj_dual_solver.hpp"
#include "inputs.hpp"

void transient_solve(int argc, char *argv[], hephaestus::Inputs inputs);
