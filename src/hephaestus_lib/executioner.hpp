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

class TransientExecutioner {
private:
  mutable double t_step; // Time step
  double t_initial;      // Start time
  double t_final;        // End time
  mutable double t;      // Current time
  int vis_steps;         // Number of cyces between each output update
  bool visualization; // Flag to control whether GLVis visualisation is required
  mutable bool last_step; // Flag to check if current step is final
  int myid_;
  int num_procs_;

  mfem::ParMesh *pmesh;
  int order;
  hephaestus::BCMap *bc_map;
  hephaestus::DomainProperties *domain_properties;
  hephaestus::AuxKernels *auxkernels;
  hephaestus::Postprocessors *postprocessors;
  hephaestus::Sources *sources;
  hephaestus::Outputs *outputs;
  std::map<std::string, mfem::DataCollection *> *data_collections;
  hephaestus::InputParameters *solver_options;

  mfem::ODESolver *ode_solver;
  mfem::BlockVector *F;

public:
  TransientExecutioner() = default;
  explicit TransientExecutioner(const hephaestus::InputParameters &params);

  void Init();

  void Step(double dt, int it) const;

  void Solve() const;

  // Enable output to GLVis
  void EnableVisualisation() { visualization = true; };
  hephaestus::FESpaces *fespaces;
  hephaestus::GridFunctions *gridfunctions;
  hephaestus::TransientFormulation *formulation;
};

} // namespace hephaestus
