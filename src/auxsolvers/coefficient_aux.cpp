#include "coefficient_aux.hpp"

#include <utility>

namespace hephaestus
{

CoefficientAux::CoefficientAux(std::string gf_name,
                               std::string coef_name,
                               hephaestus::InputParameters solver_options)
  : _gf_name(std::move(gf_name)),
    _coef_name(std::move(coef_name)),
    _solver_options(std::move(solver_options))
{
}

void
CoefficientAux::Init(const hephaestus::GridFunctions & gridfunctions,
                     hephaestus::Coefficients & coefficients)
{
  _gf = gridfunctions.Get(_gf_name);
  _coef = coefficients._scalars.Get(_coef_name);
  _test_fes = _gf->ParFESpace();

  BuildBilinearForm();
  BuildLinearForm();
  _a_mat = std::unique_ptr<mfem::HypreParMatrix>(_a->ParallelAssemble());
  _solver = std::make_unique<hephaestus::DefaultJacobiPCGSolver>(_solver_options, *_a_mat);
}

void
CoefficientAux::Update()
{
  // NB: assume _gf, _coef, _test_fes already updated earlier in update chain.

  // Update _a.
  _a->Update();
  _a->Assemble();
  _a->Finalize();

  // Update _b.
  _b->Update();
  _b->Assemble();

  // Rebuild Hypre matrix.
  _a_mat.reset();
  _a_mat = std::unique_ptr<mfem::HypreParMatrix>(_a->ParallelAssemble());

  // Rebuild solver.
  _solver.reset();
  _solver = std::make_unique<hephaestus::DefaultJacobiPCGSolver>(_solver_options, *_a_mat);
}

void
CoefficientAux::BuildBilinearForm()
{
  _a = std::make_unique<mfem::ParBilinearForm>(_test_fes);
  _a->AddDomainIntegrator(new mfem::MassIntegrator());
  _a->Assemble();
  _a->Finalize();
}

void
CoefficientAux::BuildLinearForm()
{
  _b = std::make_unique<mfem::ParLinearForm>(_test_fes);
  _b->AddDomainIntegrator(new mfem::DomainLFIntegrator(*_coef));
  _b->Assemble();
}

void
CoefficientAux::Solve(double t)
{
  mfem::Vector x(_test_fes->GetTrueVSize()); // Gridfunction true DOFs
  _gf->ProjectCoefficient(*_coef);           // Initial condition
  _gf->GetTrueDofs(x);

  // Reassemble in case coef has changed
  _b->Update();
  _b->Assemble();

  _solver->Mult(*_b, x);
  _gf->SetFromTrueDofs(x);
}

} // namespace hephaestus
