#pragma once
#include "auxkernels.hpp"
#include "inputs.hpp"
#include "postprocessors.hpp"
#include "sources.hpp"
#include "variables.hpp"
#include <fstream>
#include <iostream>
#include <memory>

namespace hephaestus {

class SteadyExecutioner {
private:
  int vis_steps;         // Number of cyces between each output update
  bool visualization; // Flag to control whether GLVis visualisation is required
  mutable bool last_step; // Flag to check if current step is final

  mfem::ParMesh *pmesh;
  hephaestus::BCMap *bc_map;
  hephaestus::DomainProperties *domain_properties;
  hephaestus::AuxKernels *auxkernels;
  hephaestus::Postprocessors *postprocessors;
  hephaestus::Sources *sources;
  hephaestus::Outputs *outputs;
  std::map<std::string, mfem::DataCollection *> *data_collections;
  hephaestus::InputParameters *solver_options;

  // Shouldn't be a TransientFormulation really
  hephaestus::TransientFormulation *formulation;
  mfem::ODESolver *ode_solver;
  mfem::BlockVector *F;

public:
  SteadyExecutioner() = default;
  explicit SteadyExecutioner(const hephaestus::InputParameters &params);

  void Init(const hephaestus::InputParameters &params);

  // // Solve transient FE problem.
  // void Solve(const hephaestus::InputParameters &params) const;
  void Solve() const;

  // Enable output to GLVis
  void EnableVisualisation() { visualization = true; };
  hephaestus::Variables *variables;
};

} // namespace hephaestus
