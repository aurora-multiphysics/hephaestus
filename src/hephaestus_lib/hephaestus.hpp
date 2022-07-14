#pragma once
#include "hephaestus_hertz.hpp"
#include "hephaestus_joule.hpp"
#include "hephaestus_transient.hpp"

void run_hephaestus(int argc, char *argv[], hephaestus::Inputs inputs) {
  if (inputs.formulation == "Joule") {
    joule_solve(argc, argv, inputs);
  } else if (inputs.formulation == "Hertz") {
    hertz_solve(argc, argv, inputs);
  } else {
    transient_solve(argc, argv, inputs);
  }
}
