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

class Executioner {
protected:
  bool visualization; // Flag to control whether GLVis visualisation is required
  int myid_;
  int num_procs_;

  mfem::ParMesh *pmesh;
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
  Executioner() = default;
  explicit Executioner(const hephaestus::InputParameters &params);

  // Initialise owned objects
  virtual void Init(){};

  // Solve the current system of equations
  virtual void Solve() const {};

  // Execute solution strategy including any timestepping
  virtual void Execute() const {};

  // Enable output to GLVis
  void EnableVisualisation() { visualization = true; };
  hephaestus::FESpaces *fespaces;
  hephaestus::GridFunctions *gridfunctions;
};

} // namespace hephaestus
