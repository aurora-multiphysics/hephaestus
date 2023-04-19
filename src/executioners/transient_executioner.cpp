#include "transient_executioner.hpp"

namespace hephaestus {

TransientExecutioner::TransientExecutioner(
    const hephaestus::InputParameters &params)
    : ExecutionerBase(params), t_step(params.GetParam<float>("TimeStep")),
      t_initial(params.GetParam<float>("StartTime")),
      t_final(params.GetParam<float>("EndTime")), t(t_initial),
      vis_steps(params.GetOptionalParam<int>("VisualisationSteps", 1)),
      last_step(false) {
  // Set Formulation
  formulation =
      params.GetParam<hephaestus::TransientFormulation *>("Formulation");
  if (formulation->equation_system == NULL) {
    formulation->CreateEquationSystem();
  }
}

void TransientExecutioner::Init() {
  fespaces->Init(*pmesh);
  gridfunctions->Init(*pmesh, *fespaces);
  formulation->RegisterMissingVariables(*pmesh, *fespaces, *gridfunctions);
  formulation->RegisterAuxKernels(*gridfunctions, *auxkernels);
  formulation->RegisterCoefficients(*domain_properties);
  formulation->equation_system->Init(*gridfunctions, *fespaces, *bc_map,
                                     *domain_properties);
  auxkernels->Init(*gridfunctions, *domain_properties);
  sources->Init(*gridfunctions, *fespaces, *bc_map, *domain_properties);
  formulation->CreateTimeDomainOperator(
      *pmesh, order, *fespaces, *gridfunctions, *bc_map, *domain_properties,
      *sources, *solver_options);

  formulation->td_operator->SetVariables();
  F = new mfem::BlockVector(
      formulation->td_operator->true_offsets); // Vector of dofs
  formulation->td_operator->Init(*F);          // Set up initial conditions

  formulation->td_operator->SetTime(t);
  ode_solver = new mfem::BackwardEulerSolver;
  ode_solver->Init(*formulation->td_operator);

  postprocessors->Init(*gridfunctions, *domain_properties);
  auxkernels->Solve(t);

  // Set up DataCollections to track fields of interest.
  for (auto const &[name, dc_] : *data_collections) {
    outputs->RegisterOutputFields(dc_, pmesh, *gridfunctions);
    // Write initial fields to disk
    outputs->WriteOutputFields(dc_, t, 0);
  }

  // Initialize GLVis visualization and send the initial condition
  // by socket to a GLVis server.
  if (visualization) {
    outputs->InitializeGLVis(myid_, *gridfunctions);
    outputs->DisplayToGLVis(*gridfunctions);
  }
}

void TransientExecutioner::Step(double dt, int it) const {
  // Check if current time step is final
  if (t + dt >= t_final - dt / 2) {
    last_step = true;
  }

  // Advance time step.
  ode_solver->Step(*F, t, dt);
  auxkernels->Solve(t);

  // Output data
  if (last_step || (it % vis_steps) == 0) {

    // Output timestep summary to console
    outputs->WriteConsoleSummary(myid_, t, it);

    // Make sure all ranks have sent their 'v' solution before initiating
    // another set of GLVis connections (one from each rank):
    MPI_Barrier(pmesh->GetComm());

    // Send output fields to GLVis for visualisation
    if (visualization) {
      outputs->DisplayToGLVis(*gridfunctions);
    }

    // Save output fields at timestep to DataCollections
    for (auto const &[name, dc_] : *data_collections) {
      outputs->WriteOutputFields(dc_, t, it);
    }
    postprocessors->Update(t);
  }
}

void TransientExecutioner::Solve() const {
  // Initialise time variables
  t = t_initial;
  last_step = false;

  for (int it = 1; !last_step; it++) {
    Step(t_step, it);
  }
}

} // namespace hephaestus