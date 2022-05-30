//            ---------------------------------------------------
//            ESolver:  Low-Frequency Magnetodynamics Simulation
//            ---------------------------------------------------
//
// This miniapp solves low frequency magnetodynamics problems using the AV
// formulation.
//
// (ν∇×A, ∇×A') + (σ(dA/dt + ∇ V), A') - (J0, A') - <(ν∇×A) × n, A'> = 0
// (σ(dA/dt + ∇ V), ∇ V') + <n.J, V'> = 0

#pragma once
#include "../common/pfem_extras.hpp"
#include "e_solver.hpp"
#include "inputs.hpp"

void e_solve(int argc, char *argv[], hephaestus::Inputs inputs) {

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // Read the serial mesh from the given mesh file on all processors. We can
  // handle triangular, quadrilateral, tetrahedral and hexahedral meshes
  // with the same code.
  mfem::Mesh *mesh = new mfem::Mesh(inputs.mesh);
  int dim = mesh->Dimension();
  int max_attr = mesh->bdr_attributes.Max();
  mesh->EnsureNCMesh(); // Required for mesh refinement
  mfem::ParMesh pmesh = mfem::ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;

  int order = inputs.order;
  hephaestus::BCMap bc_map(inputs.bc_map);
  hephaestus::DomainProperties domain_properties(inputs.domain_properties);
  hephaestus::ESolver esolver(pmesh, order, bc_map, domain_properties);

  // Allocate memory to store
  int Vsize_h1 = esolver.H1FESpace_->GetVSize();
  int Vsize_nd = esolver.HCurlFESpace_->GetVSize();

  std::cout << Vsize_h1 << " ok ";
  mfem::Array<int> true_offset(3);
  true_offset[0] = 0;
  true_offset[1] = Vsize_h1;
  true_offset[2] = true_offset[1] + Vsize_nd;

  mfem::BlockVector F(true_offset);
  std::cout << Vsize_h1;

  esolver.Init(F);
  std::cout << Vsize_h1;

  mfem::ODESolver *ode_solver = new mfem::BackwardEulerSolver;
  bool visualization = true;
  bool visit = true;

  // Initialize GLVis visualization
  if (visualization) {
    esolver.InitializeGLVis();
  }

  // Initialize VisIt visualization
  mfem::VisItDataCollection visit_dc("AV-Parallel", &pmesh);

  if (visit) {
    esolver.RegisterVisItFields(visit_dc);
  }

  // Write initial fields to disk for VisIt
  if (visit) {
    esolver.WriteVisItFields(0);
  }

  // Send the initial condition by socket to a GLVis server.
  if (visualization) {
    esolver.DisplayToGLVis();
  }

  double ti = 0.0; // initial time
  double dt = 0.5;
  double t_final = 2.5;
  int vis_steps = 1;

  double t = ti;
  esolver.SetTime(t);
  ode_solver->Init(esolver);

  bool last_step = false;
  for (int it = 1; !last_step; it++) {
    if (t + dt >= t_final - dt / 2) {
      last_step = true;
    }

    // F is the vector of dofs, t is the current time, and dt is the time step
    // to advance.
    ode_solver->Step(F, t, dt);

    if (last_step || (it % vis_steps) == 0) {
      if (myid == 0) {
        std::cout << std::fixed;
        std::cout << "step " << std::setw(6) << it << ",\tt = " << std::setw(6)
                  << std::setprecision(3) << t
                  << ",\tdot(E, J) = " << std::endl;
      }

      // Make sure all ranks have sent their 'v' solution before initiating
      // another set of GLVis connections (one from each rank):
      MPI_Barrier(pmesh.GetComm());

      if (visualization) {
        esolver.DisplayToGLVis();
      }

      if (visit) {
        esolver.WriteVisItFields(it);
      }
    }
  }
  if (myid == 0) {
    std::cout << "\nSolved" << std::endl;
  }
}
