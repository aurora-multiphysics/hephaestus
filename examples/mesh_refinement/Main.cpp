#include "mfem.hpp"
#include "hephaestus.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace hephaestus;

int
main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

  const char * mesh_file = "../../data/star.mesh";

  mfem::Mesh mesh(mesh_file, 1, 1);
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

  SteadyStateEquationSystemProblemBuilder problem_builder;

  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);
  problem_builder.SetMesh(pmesh);

  problem_builder.AddFESpace("fespace", "H1_3D_P1");
  problem_builder.AddGridFunction("gridfunction", "fespace");

  // Create useful coefficients and add.
  Coefficients coefficients;
  coefficients._scalars.Register("one", std::make_shared<mfem::ConstantCoefficient>(1.0));
  coefficients._scalars.Register("zero", std::make_shared<mfem::ConstantCoefficient>(0.0));
  problem_builder.SetCoefficients(coefficients);

  // Construct the operator so we can add kernels.
  problem_builder.ConstructOperator();

  InputParameters params_bilinear;
  params_bilinear.SetParam(std::string("CoefficientName"), std::string("one"));

  auto diffusion_kernel = std::make_shared<DiffusionKernel>(params_bilinear);
  problem_builder.AddKernel<mfem::ParBilinearForm>(std::string("gridfunction"),
                                                   std::move(diffusion_kernel));

  InputParameters params_linear;
  params_linear.SetParam(std::string("CoefficientName"), std::string("one"));

  auto linear_kernel = std::make_shared<LinearKernel>(params_linear);
  problem_builder.AddKernel<mfem::ParLinearForm>(std::string("gridfunction"),
                                                 std::move(linear_kernel));

  // All boundary attributes will be used for essential (Dirichlet) BC.
  mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
  ess_bdr = 1;

  auto bc = std::make_shared<ScalarDirichletBC>(
      "gridfunction", ess_bdr, coefficients._scalars.Get("one"));
  problem_builder.AddBoundaryCondition("gridfunction", std::move(bc));

  // Finalize construction and initialization of problem.
  problem_builder.FinalizeProblem(false);

  // Finished building; return problem.
  auto ss_eqn_system_problem = problem_builder.ReturnProblem();

  auto * problem_operator = ss_eqn_system_problem->GetOperator();
  auto * equation_system = ss_eqn_system_problem->GetEquationSystem();

  // Setup solver parameters.
  hephaestus::InputParameters params;
  params.SetParam("Tolerance", 1.0e-9);
  params.SetParam("MaxIter", 100);

  const int max_iteration = 10;

  for (int it = 0; it < max_iteration; it++)
  {
    std::cout << "Now doing iteration " << it << std::endl;

    // Form linear system AX=B; remove BCs.
    equation_system->FormLinearSystem(problem_operator->_equation_system_operator,
                                      problem_operator->_true_x,
                                      problem_operator->_true_rhs);

    // Covnert operator to hypre matrix A.
    auto hypre_matrix = problem_operator->_equation_system_operator.As<mfem::HypreParMatrix>();

    // Construct solver; solve:
    hephaestus::DefaultH1PCGSolver solver(params, *hypre_matrix);
    solver.Mult(problem_operator->_true_rhs, problem_operator->_true_x);

    // Add back in boundary conditions.
    equation_system->RecoverFEMSolution(problem_operator->_true_x,
                                        ss_eqn_system_problem->_gridfunctions);

    // pmesh->UniformRefinement();

    // Call update methods, etc, rebuild equation system, update boundary conditions.
  }

  MPI_Finalize();
  return 0;
}
