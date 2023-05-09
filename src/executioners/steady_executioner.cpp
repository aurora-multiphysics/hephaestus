#include "steady_executioner.hpp"

namespace hephaestus {

SteadyExecutioner::SteadyExecutioner(const hephaestus::InputParameters &params)
    : ExecutionerBase(params) {
  // Set Formulation
  formulation = params.GetParam<hephaestus::SteadyFormulation *>("Formulation");
  // if (formulation->equation_system == NULL) {
  //   formulation->CreateEquationSystem();
  // }
}

// fd_operator
// SetVariables()
// Init(*F)
// Solve(*F)
void SteadyExecutioner::Init() {
  fespaces->Init(*pmesh);
  gridfunctions->Init(*pmesh, *fespaces);
  formulation->RegisterMissingVariables(*pmesh, *fespaces, *gridfunctions);
  formulation->RegisterAuxKernels(*gridfunctions, *auxkernels);
  formulation->RegisterCoefficients(*domain_properties);
  auxkernels->Init(*gridfunctions, *domain_properties);
  sources->Init(*gridfunctions, *fespaces, *bc_map, *domain_properties);
  formulation->CreateFrequencyDomainOperator(*pmesh, *fespaces, *gridfunctions,
                                             *bc_map, *domain_properties,
                                             *sources, *solver_options);

  formulation->fd_operator->SetVariables();
  F = new mfem::BlockVector(
      formulation->fd_operator->true_offsets); // Vector of dofs
  formulation->fd_operator->Init(*F);          // Set up initial conditions

  // ode_solver->Init(*formulation->td_operator);

  postprocessors->Init(*gridfunctions, *domain_properties);
  auxkernels->Solve();

  // Set up DataCollections to track fields of interest.
  for (auto const &[name, dc_] : *data_collections) {
    outputs->RegisterOutputFields(dc_, pmesh, *gridfunctions);
    // Write initial fields to disk
    outputs->WriteOutputFields(dc_, 0.0, 0);
  }

  // Initialize GLVis visualization and send the initial condition
  // by socket to a GLVis server.
  if (visualization) {
    outputs->InitializeGLVis(myid_, *gridfunctions);
    outputs->DisplayToGLVis(*gridfunctions);
  }
}

void SteadyExecutioner::Solve() const {
  // Advance time step.
  formulation->fd_operator->Solve(*F);
  auxkernels->Solve();

  // Output data
  // Output timestep summary to console
  outputs->WriteConsoleSummary(myid_, 1.0, 1);

  // Make sure all ranks have sent their 'v' solution before initiating
  // another set of GLVis connections (one from each rank):
  MPI_Barrier(pmesh->GetComm());

  // Send output fields to GLVis for visualisation
  if (visualization) {
    outputs->DisplayToGLVis(*gridfunctions);
  }

  // Save output fields at timestep to DataCollections
  for (auto const &[name, dc_] : *data_collections) {
    outputs->WriteOutputFields(dc_, 1.0, 1);
  }
  postprocessors->Update();
}

} // namespace hephaestus
