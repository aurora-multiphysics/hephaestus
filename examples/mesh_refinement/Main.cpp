#include "mfem.hpp"
#include "hephaestus.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int
main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  const char * mesh_file = "../../data/star.mesh";

  Mesh mesh(mesh_file, 1, 1);
  int dim = mesh.Dimension();
  int sdim = mesh.SpaceDimension();

  if (mesh.NURBSext)
  {
    for (int i = 0; i < 2; i++)
    {
      mesh.UniformRefinement();
    }
    mesh.SetCurvature(2);
  }

  hephaestus::SteadyStateProblemBuilder problem_builder;

  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);
  problem_builder.SetMesh(pmesh);

  problem_builder.AddFESpace("fespace", "H1_3D_P1");

  problem_builder.AddGridFunction("gridfunction", "fespace");

  hephaestus::Coefficients coefficients;

  coefficients._scalars.Register("one", std::make_shared<mfem::ConstantCoefficient>(1.0));
  coefficients._scalars.Register("zero", std::make_shared<mfem::ConstantCoefficient>(0.0));
  problem_builder.SetCoefficients(coefficients);

  problem_builder.ConstructOperator();
  problem_builder.ConstructEquationSystem();

  hephaestus::InputParameters params_bilinear;
  params_bilinear.SetParam(std::string("CoefficientName"), std::string("one"));

  auto diffusion_kernel = std::make_shared<hephaestus::DiffusionKernel>(params_bilinear);
  problem_builder.AddKernel<mfem::ParBilinearForm>(std::string("gridfunction"),
                                                   std::move(diffusion_kernel));

  hephaestus::InputParameters params_linear;
  params_linear.SetParam(std::string("CoefficientName"), std::string("one"));

  auto linear_kernel = std::make_shared<hephaestus::LinearKernel>(params_linear);
  problem_builder.AddKernel<mfem::ParLinearForm>(std::string("gridfunction"),
                                                 std::move(linear_kernel));

  // 8. All boundary attributes will be used for essential (Dirichlet) BC.
  mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
  ess_bdr = 1;

  auto bc = std::make_shared<hephaestus::ScalarDirichletBC>(
      "gridfunction", ess_bdr, coefficients._scalars.Get("one"));
  problem_builder.AddBoundaryCondition("gridfunction", std::move(bc));

  // Register everything.
  problem_builder.RegisterFESpaces();
  problem_builder.RegisterGridFunctions();
  problem_builder.RegisterAuxSolvers();
  problem_builder.RegisterCoefficients();

  // Initialize.
  problem_builder.ConstructState();
  problem_builder.InitializeAuxSolvers();
  problem_builder.InitializeKernels();
  problem_builder.InitializeOutputs();

  auto problem = problem_builder.ReturnProblem();

  problem->GetOperator()->SetGridFunctions();

  // Construct eqn system.
  problem->GetEquationSystem()->Init(
      problem->_gridfunctions,
      problem->_fespaces,
      problem->_bc_map,
      problem->_coefficients); // Becuase not initialized properly for steady state problems.
  problem->GetEquationSystem()->BuildEquationSystem(problem->_bc_map, problem->_sources);

  const int max_it = 10;

  for (int it = 0; it < max_it; it++)
  {
    std::cout << "Now doing iteration " << it << std::endl;

    problem->GetEquationSystem()->FormLinearSystem(
        problem->GetOperator()->_equation_system_operator,
        problem->GetOperator()->_true_x,
        problem->GetOperator()->_true_rhs);

    // // Solve.
    // GSSmoother M((SparseMatrix &)(*A));
    // PCG(*A, M, B, X, 3, 200, 1e-12, 0.0);

    // problem->GetEquationSystem()->RecoverFEMSolution(X, problem->_gridfunctions);

    // pmesh->UniformRefinement();

    // Call update methods, etc, rebuild equation system, update boundary conditions.
  }

  MPI_Finalize();
  return 0;
}
