#include "transient_executioner.hpp"

namespace hephaestus {

TransientExecutioner::TransientExecutioner(
    const hephaestus::InputParameters &params)
    : Executioner(params), t_step(params.GetParam<float>("TimeStep")),
      t_initial(params.GetParam<float>("StartTime")),
      t_final(params.GetParam<float>("EndTime")), t(t_initial), it(0),
      vis_steps(params.GetOptionalParam<int>("VisualisationSteps", 1)),
      last_step(false) {
  // Set Formulation
  formulation =
      params.GetParam<hephaestus::TransientFormulation *>("Formulation");
  td_equation_system = dynamic_cast<hephaestus::TimeDependentEquationSystem *>(
      formulation->CreateEquationSystem());
}

void TransientExecutioner::Init() {
  fespaces->Init(*pmesh);
  gridfunctions->Init(*pmesh, *fespaces);
  formulation->RegisterMissingVariables(*pmesh, *fespaces, *gridfunctions);
  formulation->RegisterAuxKernels(*gridfunctions, *auxkernels);
  formulation->RegisterCoefficients(*domain_properties);
  td_equation_system->Init(*gridfunctions, *fespaces, *bc_map,
                           *domain_properties);
  auxkernels->Init(*gridfunctions, *domain_properties);
  sources->Init(*gridfunctions, *fespaces, *bc_map, *domain_properties);
  td_operator = dynamic_cast<hephaestus::TimeDomainEquationSystemOperator *>(
      formulation->CreateOperator(*pmesh, *fespaces, *gridfunctions, *bc_map,
                                  *domain_properties, *sources,
                                  *solver_options));
  td_operator->SetEquationSystem(td_equation_system);

  td_operator->SetVariables();
  F = new mfem::BlockVector(td_operator->true_offsets); // Vector of dofs
  td_operator->Init(*F); // Set up initial conditions

  td_operator->SetTime(t);
  ode_solver = new mfem::BackwardEulerSolver;
  ode_solver->Init(*td_operator);

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

void TransientExecutioner::Solve() const { Step(t_step, it); }

void TransientExecutioner::Execute() const {
  // Initialise time variables
  t = t_initial;
  last_step = false;

  for (it = 1; !last_step; it++) {
    Solve();
  }
}

} // namespace hephaestus
