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

  Coefficients coefficients;
  coefficients._scalars.Register("one", std::make_shared<mfem::ConstantCoefficient>(1.0));
  coefficients._scalars.Register("zero", std::make_shared<mfem::ConstantCoefficient>(0.0));

  problem_builder.SetCoefficients(coefficients);

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

  // 8. All boundary attributes will be used for essential (Dirichlet) BC.
  mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
  ess_bdr = 1;

  auto bc = std::make_shared<ScalarDirichletBC>(
      "gridfunction", ess_bdr, coefficients._scalars.Get("one"));
  problem_builder.AddBoundaryCondition("gridfunction", std::move(bc));

  problem_builder.FinalizeProblem(false);

  auto steady_state_problem = problem_builder.ReturnProblem();

  auto steady_state_equation_system_problem = std::unique_ptr<SteadyStateEquationSystemProblem>(
      static_cast<SteadyStateEquationSystemProblem *>(steady_state_problem.release()));

  EquationSystemProblemOperator * problem_operator =
      steady_state_equation_system_problem->GetOperator();
  EquationSystem * equation_system = steady_state_equation_system_problem->GetEquationSystem();

  const int max_it = 10;

  for (int it = 0; it < max_it; it++)
  {
    std::cout << "Now doing iteration " << it << std::endl;

    equation_system->FormLinearSystem(problem_operator->_equation_system_operator,
                                      problem_operator->_true_x,
                                      problem_operator->_true_rhs);

    // Solve.
    mfem::GSSmoother smoother((mfem::SparseMatrix &)(*problem_operator->_equation_system_operator));
    mfem::PCG(*problem_operator->_equation_system_operator,
              smoother,
              problem_operator->_true_rhs,
              problem_operator->_true_x,
              3,
              200,
              1e-12,
              0.0);

    equation_system->RecoverFEMSolution(problem_operator->_true_x,
                                        steady_state_equation_system_problem->_gridfunctions);

    // pmesh->UniformRefinement();

    // Call update methods, etc, rebuild equation system, update boundary conditions.
  }

  MPI_Finalize();
  return 0;
}
