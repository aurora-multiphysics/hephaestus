// Solves a time-domain formulation->
// TODO: Add to Executioners?

#pragma once
#include "../common/pfem_extras.hpp"
#include "factory.hpp"
#include "inputs.hpp"

void transient_solve(int argc, char *argv[], hephaestus::Inputs inputs);
