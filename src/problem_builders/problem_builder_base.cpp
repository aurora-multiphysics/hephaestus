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
  GetProblem()->_solvers.SetComm(pmesh->GetComm());
  MPI_Comm_size(pmesh->GetComm(), &(GetProblem()->_num_procs));
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
  GetProblem()->_solvers.SetSolverOptions(solver_options);
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

void
ProblemSolvers::SetSolver(std::shared_ptr<mfem::Solver> & solver,
                          std::string solve_type,
                          mfem::Operator & A)
{

  if (!_comm_set)
  {
    mfem::mfem_error("MPI comm must be set in the solvers class before setting up solver!");
  }

  if (solve_type == "PCG")
  {
    std::shared_ptr<mfem::HyprePCG> solver_buffer{std::make_shared<mfem::HyprePCG>(_comm)};

    solver_buffer->SetTol(_solver_options.GetParam<float>("Tolerance"));
    solver_buffer->SetAbsTol(_solver_options.GetParam<float>("AbsTolerance"));
    solver_buffer->SetMaxIter(_solver_options.GetParam<unsigned int>("MaxIter"));
    solver_buffer->SetPrintLevel(_solver_options.GetParam<int>("PrintLevel"));
    solver_buffer->SetOperator(A);

    if (_linear_preconditioner)
      solver_buffer->SetPreconditioner(
          *std::dynamic_pointer_cast<mfem::HypreSolver>(_linear_preconditioner));

    solver = solver_buffer;
  }
  else if (solve_type == "GMRES")
  {
    std::shared_ptr<mfem::HypreGMRES> solver_buffer{std::make_shared<mfem::HypreGMRES>(_comm)};

    solver_buffer->SetTol(_solver_options.GetParam<float>("Tolerance"));
    solver_buffer->SetAbsTol(_solver_options.GetParam<float>("AbsTolerance"));
    solver_buffer->SetMaxIter(_solver_options.GetParam<unsigned int>("MaxIter"));
    solver_buffer->SetPrintLevel(_solver_options.GetParam<int>("PrintLevel"));
    solver_buffer->SetOperator(A);

    if (_linear_preconditioner)
      solver_buffer->SetPreconditioner(
          *std::dynamic_pointer_cast<mfem::HypreSolver>(_linear_preconditioner));

    solver = solver_buffer;
  }
  else if (solve_type == "FGMRES")
  {
    std::shared_ptr<mfem::HypreFGMRES> solver_buffer{std::make_shared<mfem::HypreFGMRES>(_comm)};

    solver_buffer->SetTol(_solver_options.GetParam<float>("Tolerance"));
    solver_buffer->SetMaxIter(_solver_options.GetParam<unsigned int>("MaxIter"));
    solver_buffer->SetPrintLevel(_solver_options.GetParam<int>("PrintLevel"));
    solver_buffer->SetOperator(A);

    if (_linear_preconditioner)
      solver_buffer->SetPreconditioner(
          *std::dynamic_pointer_cast<mfem::HypreSolver>(_linear_preconditioner));

    solver = solver_buffer;
  }
  else if (solve_type == "AMS")
  {
    if (!_edge_fespace)
      mfem::mfem_error("_edge_fespace must be set within ProblemSolvers class before calling AMS!");

    std::shared_ptr<mfem::HypreAMS> solver_buffer{std::make_shared<mfem::HypreAMS>(_edge_fespace)};

    solver_buffer->SetSingularProblem();
    solver_buffer->SetPrintLevel(_solver_options.GetParam<int>("PrintLevel"));
    solver_buffer->SetOperator(A);

    solver = solver_buffer;
  }
  else if (solve_type == "AMG")
  {
    std::shared_ptr<mfem::HypreBoomerAMG> solver_buffer{std::make_shared<mfem::HypreBoomerAMG>()};

    // solver_buffer->SetTol(_solver_options.GetParam<float>("Tolerance"));
    // solver_buffer->SetMaxIter(_solver_options.GetParam<unsigned int>("MaxIter"));
    solver_buffer->SetPrintLevel(_solver_options.GetParam<int>("PrintLevel"));
    solver_buffer->SetOperator(A);

    solver = solver_buffer;
  }
  else if (solve_type == "SuperLU")
  {
    std::shared_ptr<mfem::SuperLUSolver> solver_buffer{
        std::make_shared<mfem::SuperLUSolver>(MPI_COMM_WORLD)};

    mfem::SuperLURowLocMatrix a_super_lu(A);
    solver_buffer->SetOperator(a_super_lu);

    solver = solver_buffer;
  }
  else
  {
    mfem::mfem_error(("Solver or Preconditioner type " + solve_type +
                      " not recognised. Please set to one of the following options: PCG, GMRES, "
                      "FGMRES, AMS, AMG, SuperLU.")
                         .c_str());
  }

  // Add others!
}

void
ProblemSolvers::SetLinearPreconditioner(mfem::Operator & A)
{
  if (!_solver_options.HasParam("LinearPreconditioner"))
  {
    std::cout
        << "Warning: LinearPreconditioner field not set. Proceeding without a preconditioner\n";
  }
  else
  {
    auto solve_type = _solver_options.GetParam<std::string>("LinearPreconditioner");
    SetSolver(_linear_preconditioner, solve_type, A);
  }
}

void
ProblemSolvers::SetLinearSolver(mfem::Operator & A)
{

  if (!_solver_options.HasParam("LinearSolver"))
  {
    mfem::mfem_error("LinearSolver field not set. Please set to one of the following options: PCG, "
                     "GMRES, FGMRES, AMS, AMG, SuperLU.");
  }
  else
  {
    auto solve_type = _solver_options.GetParam<std::string>("LinearSolver");
    SetSolver(_linear_solver, solve_type, A);
  }
}

void
ProblemSolvers::SetComm(const MPI_Comm comm)
{
  _comm = comm;
  _comm_set = true;
}

void
ProblemSolvers::SetSolverOptions(const hephaestus::InputParameters solver_options)
{
  _solver_options = solver_options;
}

void
ProblemSolvers::SetEdgeFESpace(mfem::ParFiniteElementSpace * edge_fes)
{
  _edge_fespace = edge_fes;
}

} // namespace hephaestus
