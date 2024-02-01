#include "problem_builder.hpp"

namespace hephaestus
{

Problem::~Problem()
{
  // Ensure that all owned memory is properly freed!
  _f.reset();
  _ode_solver.reset();
}

void
ProblemBuilder::SetMesh(std::shared_ptr<mfem::ParMesh> pmesh)
{
  GetProblem()->_pmesh = pmesh;
  MPI_Comm_rank(pmesh->GetComm(), &(GetProblem()->_myid));
}

void
ProblemBuilder::SetFESpaces(hephaestus::FESpaces & fespaces)
{
  GetProblem()->_fespaces = fespaces;
}

void
ProblemBuilder::SetGridFunctions(hephaestus::GridFunctions & gridfunctions)
{
  GetProblem()->_gridfunctions = gridfunctions;
}

void
ProblemBuilder::SetBoundaryConditions(hephaestus::BCMap & bc_map)
{
  GetProblem()->_bc_map = bc_map;
}

void
ProblemBuilder::SetAuxSolvers(hephaestus::AuxSolvers & preprocessors)
{
  GetProblem()->_preprocessors = preprocessors;
}

void
ProblemBuilder::SetPostprocessors(hephaestus::AuxSolvers & postprocessors)
{
  GetProblem()->_postprocessors = postprocessors;
}

void
ProblemBuilder::SetSources(hephaestus::Sources & sources)
{
  GetProblem()->_sources = sources;
}

void
ProblemBuilder::SetOutputs(hephaestus::Outputs & outputs)
{
  GetProblem()->_outputs = outputs;
}

void
ProblemBuilder::SetSolverOptions(hephaestus::InputParameters & solver_options)
{
  GetProblem()->_solver_options = solver_options;
}

void
ProblemBuilder::SetSolverOperator(std::unique_ptr<mfem::Solver> & solver,
                                  std::string solve_type,
                                  mfem::HypreParMatrix & A)
{

  if (solve_type == "PCG")
  {
    solver = std::make_unique<mfem::HyprePCG>(A);
  }
  else if (solve_type == "GMRES")
  {
    solver = std::make_unique<mfem::HypreGMRES>(A);
    // GetProblem()->_jacobian_preconditioner->SetKDim(GetProblem()->_solver_options.GetParam<int>("KDim"));
  }
  else if (solve_type == "FGMRES")
  {
    solver = std::make_unique<mfem::HypreFGMRES>(A);
    // GetProblem()->_jacobian_preconditioner->SetKDim(GetProblem()->_solver_options.GetParam<int>("KDim"));
  }
  // Add others!
}

void
ProblemBuilder::SetJacobianPreconditioner(mfem::HypreParMatrix & A)
{

  auto solve_type = GetProblem()->_solver_options.GetParam<std::string>("JacobianPreconditioner");

  SetSolverOperator(GetProblem()->_jacobian_preconditioner, solve_type, A);

  if (!GetProblem()->_jacobian_preconditioner)
    std::cout
        << "Warning: Jacobian preconditioner type not recognised or not set. Proceeding without "
           "preconditioner. If you wish to use a preconditioner, please set "
           "JacobianPreconditioner field to one of the following values: PCG, GMRES, FGMRES"
        << std::endl;

  // if (GetProblem()->_jacobian_preconditioner)
  //{
  //   GetProblem()->_jacobian_preconditioner->SetTol(GetProblem()->_solver_options.GetParam<float>("Tolerance"));
  //   GetProblem()->_jacobian_preconditioner->SetAbsTol(GetProblem()->_solver_options.GetParam<float>("AbsTolerance"));
  //   GetProblem()->_jacobian_preconditioner->SetMaxIter(GetProblem()->_solver_options.GetParam<int>("MaxIter"));
  //   GetProblem()->_jacobian_preconditioner->SetPrintLevel(GetProblem()->_solver_options.GetParam<int>("PrintLevel"));
  // }
}

void
ProblemBuilder::SetJacobianSolver(mfem::HypreParMatrix & A)
{

  auto solve_type = GetProblem()->_solver_options.GetParam<std::string>("JacobianSolver");

  SetSolverOperator(GetProblem()->_jacobian_solver, solve_type, A);

  if (!GetProblem()->_jacobian_solver)
    mfem::mfem_error("Jacobian solver type not recognised! Please set JacobianSolver field "
                     "to one of the following values: PCG, GMRES, FGMRES");

  // GetProblem()->_jacobian_solver->SetTol(GetProblem()->_solver_options.GetParam<float>("Tolerance"));
  // GetProblem()->_jacobian_solver->SetAbsTol(GetProblem()->_solver_options.GetParam<float>("AbsTolerance"));
  // GetProblem()->_jacobian_solver->SetMaxIter(GetProblem()->_solver_options.GetParam<int>("MaxIter"));
  // GetProblem()->_jacobian_solver->SetPrintLevel(GetProblem()->_solver_options.GetParam<int>("PrintLevel"));
  // if (GetProblem()->_jacobian_preconditioner)
  //   GetProblem()->_jacobian_solver->SetPreconditioner(GetProblem()->_jacobian_preconditioner);
}

void
ProblemBuilder::SetFEMPreconditioner(mfem::HypreParMatrix & A)
{

  auto solve_type = GetProblem()->_solver_options.GetParam<std::string>("FEMPreconditioner");

  SetSolverOperator(GetProblem()->_fem_preconditioner, solve_type, A);

  if (!GetProblem()->_fem_preconditioner)
    std::cout << "Warning: FEM preconditioner type not recognised or not set. Proceeding without "
                 "preconditioner. If you wish to use a preconditioner, please set "
                 "FEMPreconditioner field to one of the following values: PCG, GMRES, FGMRES"
              << std::endl;

  // if (GetProblem()->_fem_preconditioner)
  //{
  //   GetProblem()->_fem_preconditioner->SetTol(GetProblem()->_solver_options.GetParam<float>("Tolerance"));
  //   GetProblem()->_fem_preconditioner->SetAbsTol(GetProblem()->_solver_options.GetParam<float>("AbsTolerance"));
  //   GetProblem()->_fem_preconditioner->SetMaxIter(GetProblem()->_solver_options.GetParam<int>("MaxIter"));
  //   GetProblem()->_fem_preconditioner->SetPrintLevel(GetProblem()->_solver_options.GetParam<int>("PrintLevel"));
  // }
}

void
ProblemBuilder::SetFEMSolver(mfem::HypreParMatrix & A)
{

  auto solve_type = GetProblem()->_solver_options.GetParam<std::string>("FEMSolver");

  SetSolverOperator(GetProblem()->_fem_solver, solve_type, A);

  if (!GetProblem()->_fem_solver)
    mfem::mfem_error("FEM solver type not recognised! Please set FEMSolver field "
                     "to one of the following values: PCG, GMRES, FGMRES");

  // GetProblem()->_fem_solver->SetTol(GetProblem()->_solver_options.GetParam<float>("Tolerance"));
  // GetProblem()->_fem_solver->SetAbsTol(GetProblem()->_solver_options.GetParam<float>("AbsTolerance"));
  // GetProblem()->_fem_solver->SetMaxIter(GetProblem()->_solver_options.GetParam<int>("MaxIter"));
  // GetProblem()->_fem_solver->SetPrintLevel(GetProblem()->_solver_options.GetParam<int>("PrintLevel"));
  // if (GetProblem()->_fem_preconditioner)
  //   GetProblem()->_fem_solver->SetPreconditioner(GetProblem()->_fem_preconditioner);
}

void
ProblemBuilder::SetCoefficients(hephaestus::Coefficients & coefficients)
{
  GetProblem()->_coefficients = coefficients;
}

void
ProblemBuilder::AddFESpace(std::string fespace_name, std::string fec_name, int vdim, int ordering)
{
  if (GetProblem()->_fespaces.Has(fespace_name))
  {
    const std::string error_message = "A fespace with the name " + fespace_name +
                                      " has already been added to the problem fespaces.";
    mfem::mfem_error(error_message.c_str());
  }
  if (!GetProblem()->_fecs.Has(fec_name))
  {
    mfem::FiniteElementCollection * fec = mfem::FiniteElementCollection::New(fec_name.c_str());
    GetProblem()->_fecs.Register(fec_name, fec, true);
  }

  if (!GetProblem()->_fespaces.Has(fespace_name))
  {
    mfem::ParMesh * pmesh = GetProblem()->_pmesh.get();
    if (pmesh == nullptr)
    {
      MFEM_ABORT("ParMesh not found when trying to add " << fespace_name << " to fespaces.");
    }
    auto pfes = std::make_shared<mfem::ParFiniteElementSpace>(
        GetProblem()->_pmesh.get(), GetProblem()->_fecs.Get(fec_name), vdim, ordering);

    GetProblem()->_fespaces.Register(fespace_name, std::move(pfes));
  }
}

void
ProblemBuilder::AddGridFunction(std::string gridfunction_name, std::string fespace_name)
{
  if (GetProblem()->_gridfunctions.Has(gridfunction_name))
  {
    const std::string error_message = "A gridfunction with the name " + gridfunction_name +
                                      " has already been added to the problem gridfunctions.";
    mfem::mfem_error(error_message.c_str());
  }
  mfem::ParFiniteElementSpace * fespace(GetProblem()->_fespaces.Get(fespace_name));
  if (fespace == nullptr)
  {
    MFEM_ABORT("FESpace " << fespace_name << " not found in fespaces when trying to add "
                          << gridfunction_name
                          << " associated with it into gridfunctions. Please add " << fespace_name
                          << " to fespaces before adding this gridfunction.");
  }
  auto gridfunc = std::make_shared<mfem::ParGridFunction>(fespace);
  *gridfunc = 0.0;

  GetProblem()->_gridfunctions.Register(gridfunction_name, std::move(gridfunc));
}

void
ProblemBuilder::AddBoundaryCondition(std::string bc_name,
                                     std::shared_ptr<hephaestus::BoundaryCondition> bc)
{
  if (GetProblem()->_bc_map.Has(bc_name))
  {
    const std::string error_message = "A boundary condition with the name " + bc_name +
                                      " has already been added to the problem boundary conditions.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->_bc_map.Register(bc_name, std::move(bc));
}

void
ProblemBuilder::AddAuxSolver(std::string auxsolver_name, hephaestus::AuxSolver * aux, bool own_data)
{
  if (GetProblem()->_preprocessors.Has(auxsolver_name))
  {
    const std::string error_message = "An auxsolver with the name " + auxsolver_name +
                                      " has already been added to the problem preprocessors.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->_preprocessors.Register(auxsolver_name, aux, own_data);
}

void
ProblemBuilder::AddAuxSolver(std::string auxsolver_name, std::shared_ptr<hephaestus::AuxSolver> aux)
{
  if (GetProblem()->_preprocessors.Has(auxsolver_name))
  {
    const std::string error_message = "An auxsolver with the name " + auxsolver_name +
                                      " has already been added to the problem preprocessors.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->_preprocessors.Register(auxsolver_name, std::move(aux));
}

void
ProblemBuilder::AddPostprocessor(std::string auxsolver_name,
                                 std::shared_ptr<hephaestus::AuxSolver> aux)
{
  if (GetProblem()->_postprocessors.Has(auxsolver_name))
  {
    const std::string error_message = "An auxsolver with the name " + auxsolver_name +
                                      " has already been added to the problem postprocessors.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->_postprocessors.Register(auxsolver_name, std::move(aux));
}

void
ProblemBuilder::AddSource(std::string source_name, std::shared_ptr<hephaestus::Source> source)
{
  if (GetProblem()->_sources.Has(source_name))
  {
    const std::string error_message =
        "A source with the name " + source_name + " has already been added to the problem sources.";
    mfem::mfem_error(error_message.c_str());
  }
  GetProblem()->_sources.Register(source_name, std::move(source));
}

void
ProblemBuilder::InitializeAuxSolvers()
{
  GetProblem()->_preprocessors.Init(GetProblem()->_gridfunctions, GetProblem()->_coefficients);
  GetProblem()->_postprocessors.Init(GetProblem()->_gridfunctions, GetProblem()->_coefficients);
}

void
ProblemBuilder::InitializeOutputs()
{
  GetProblem()->_outputs.Init(GetProblem()->_gridfunctions);
}

} // namespace hephaestus
