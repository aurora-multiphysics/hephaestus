#include "helmholtz_projector.hpp"

namespace hephaestus {

HelmholtzProjector::HelmholtzProjector(
    const hephaestus::InputParameters &params)
    : h1_fespace_name_(
          params.GetOptionalParam<std::string>("H1FESpaceName", "H1FES_Name")),
      hcurl_fespace_name_(params.GetOptionalParam<std::string>(
          "HCurlFESpaceName", "HCurlFES_Name")),
      gf_grad_name_(params.GetParam<std::string>("VectorGridFunctionName")),
      gf_name_(params.GetOptionalParam<std::string>("ScalarGridFunctionName",
                                                    "ScalarGF_Name")),
      H1FESpace_(nullptr), HCurlFESpace_(nullptr), q_(nullptr), g(nullptr),
      div_free_src_gf_(nullptr), gDiv_(nullptr), weakDiv_(nullptr),
      grad_(nullptr), a0_(nullptr) {

  hephaestus::InputParameters default_pars;
  default_pars.SetParam("Tolerance", float(1.0e-20));
  default_pars.SetParam("AbsTolerance", float(1.0e-20));
  default_pars.SetParam("MaxIter", (unsigned int)1000);
  default_pars.SetParam("PrintLevel", 1);

  solver_options_ = params.GetOptionalParam<hephaestus::InputParameters>(
        "SolverOptions", default_pars);

  }

HelmholtzProjector::~HelmholtzProjector() {

  if (gDiv_ != nullptr)
    delete gDiv_;
  if (weakDiv_ != nullptr)
    delete weakDiv_;
  if (grad_ != nullptr)
    delete grad_;
  if (a0_ != nullptr)
    delete a0_;
}

void HelmholtzProjector::Project(hephaestus::GridFunctions &gridfunctions,
                                 const hephaestus::FESpaces &fespaces,
                                 hephaestus::BCMap &bc_map) {

  // Retrieving vector GridFunction. This is the only mandatory one
  div_free_src_gf_ = gridfunctions.Get(gf_grad_name_);
  if (div_free_src_gf_ == nullptr) {
    const std::string error_message = gf_grad_name_ +
                                      " not found in gridfunctions when "
                                      "creating HelmholtzProjector\n";
    mfem::mfem_error(error_message.c_str());
  }

  HCurlFESpace_ = fespaces.Get(hcurl_fespace_name_);
  if (HCurlFESpace_ == nullptr) {
    std::cout << hcurl_fespace_name_ + " not found in fespaces when "
                                       "creating HelmholtzProjector. "
                                       "Obtaining from vector GridFunction.\n";
    HCurlFESpace_ = div_free_src_gf_->ParFESpace();
  }

  H1FESpace_ = fespaces.Get(h1_fespace_name_);
  if (H1FESpace_ == nullptr) {
    std::cout << h1_fespace_name_ + " not found in fespaces when "
                                    "creating HelmholtzProjector. "
                                    " Extracting from GridFunction\n";

    // Creates an H1 FES on the same mesh and with the same order as the HCurl
    // FES
    H1FESpace_ = new mfem::ParFiniteElementSpace(
        HCurlFESpace_->GetParMesh(),
        new mfem::H1_FECollection(HCurlFESpace_->GetMaxElementOrder(),
                                  HCurlFESpace_->GetParMesh()->Dimension()));
  }

  q_ = gridfunctions.Get(gf_name_);
  if (q_ == nullptr) {
    std::cout << gf_name_ + " not found in gridfunctions when "
                            "creating HelmholtzProjector. "
                            "Creating new GridFunction\n";
    q_ = new mfem::ParGridFunction(H1FESpace_);
  }

  g = new mfem::ParGridFunction(HCurlFESpace_);
  *g = *div_free_src_gf_;
  *q_ = 0.0;

  bc_map_ = &bc_map;

  setForms();
  setGrad();
  setBCs();
  solveLinearSystem();

  // Compute the irrotational component of g
  // P(g) = g - ∇Q
  grad_->Mult(*q_, *div_free_src_gf_);
  *div_free_src_gf_ -= *g;
  *div_free_src_gf_ *= -1.0;

  delete g;
  if (!gridfunctions.Has(gf_name_))
    delete q_;
  if (!fespaces.Has(h1_fespace_name_))
    delete H1FESpace_;
}

void HelmholtzProjector::setForms() {

  if (gDiv_ == nullptr)
    gDiv_ = new mfem::ParLinearForm(H1FESpace_);

  if (weakDiv_ == nullptr) {
    weakDiv_ = new mfem::ParMixedBilinearForm(HCurlFESpace_, H1FESpace_);
    weakDiv_->AddDomainIntegrator(new mfem::VectorFEWeakDivergenceIntegrator);
    weakDiv_->Assemble();
    weakDiv_->Finalize();
  }

  if (a0_ == nullptr) {
    a0_ = new mfem::ParBilinearForm(H1FESpace_);
    a0_->AddDomainIntegrator(new mfem::DiffusionIntegrator);
    a0_->Assemble();
    a0_->Finalize();
  }
}

void HelmholtzProjector::setGrad() {

  if (grad_ == nullptr) {
    grad_ = new mfem::ParDiscreteLinearOperator(H1FESpace_, HCurlFESpace_);
    grad_->AddDomainInterpolator(new mfem::GradientInterpolator());
    grad_->Assemble();
    grad_->Finalize();
  }
}

void HelmholtzProjector::setBCs() {

  // Begin Divergence-free projection
  // (g, ∇q) - (∇Q, ∇q) - <P(g).n, q> = 0
  int myid = H1FESpace_->GetMyRank();

  // <P(g).n, q>
  bc_map_->applyEssentialBCs(gf_name_, ess_bdr_tdofs_, *q_,
                             (H1FESpace_->GetParMesh()));
  bc_map_->applyIntegratedBCs(gf_name_, *gDiv_, (H1FESpace_->GetParMesh()));

  // Apply essential BC. Necessary to ensure potential at least one point is
  // fixed.
  int localsize = ess_bdr_tdofs_.Size();
  int fullsize;
  MPI_Allreduce(&fullsize, &localsize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (fullsize == 0 && myid == 0) {
    ess_bdr_tdofs_.SetSize(1);
    ess_bdr_tdofs_[0] = 0;
  }
}

void HelmholtzProjector::solveLinearSystem() {

  gDiv_->Assemble();

  // Compute the divergence of g
  // (g, ∇q)
  weakDiv_->AddMult(*g, *gDiv_, -1.0);

  // Form linear system
  // (g, ∇q) - (∇Q, ∇q) - <P(g).n, q> = 0
  // (∇Q, ∇q) = (g, ∇q) - <P(g).n, q>
  mfem::HypreParMatrix A0;
  mfem::Vector X0;
  mfem::Vector B0;
  a0_->FormLinearSystem(ess_bdr_tdofs_, *q_, *gDiv_, A0, X0, B0);

  hephaestus::DefaultGMRESSolver a0_solver(solver_options_, A0);

  a0_solver.Mult(B0, X0);
  a0_->RecoverFEMSolution(X0, *gDiv_, *q_);
}

} // namespace hephaestus