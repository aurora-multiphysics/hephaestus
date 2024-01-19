#pragma once
#include "helmholtz_projector.hpp"
#include "scalar_potential_source.hpp"
#include "source_base.hpp"

namespace hephaestus
{

double highV(const mfem::Vector & x, double t);
double lowV(const mfem::Vector & x, double t);

double calcFlux(mfem::GridFunction * v_field, int face_attr);

double calcFlux(mfem::GridFunction * v_field, int face_attr, mfem::Coefficient & q);

void inheritBdrAttributes(const mfem::ParMesh * parent_mesh, mfem::ParSubMesh * child_mesh);

// Applies the HelmholtzProjector onto the J GridFunction to clean it of any
// divergences. This is for the simplest case with no BCs
void cleanDivergence(mfem::ParGridFunction & Vec_GF, hephaestus::InputParameters solve_pars);

// The more complicated case where BCs are needed
void cleanDivergence(const hephaestus::GridFunctions & gfs,
                     const hephaestus::BCMap & bcs,
                     const std::string vec_gf_name,
                     const std::string scalar_gf_name,
                     hephaestus::InputParameters solve_pars);

void attrToMarker(const mfem::Array<int> attr_list, mfem::Array<int> & marker_list, int max_attr);

class OpenCoilSolver : public hephaestus::Source
{

public:
  OpenCoilSolver(const hephaestus::InputParameters & params,
                 mfem::Array<int>  coil_dom,
                 const std::pair<int, int> electrodes);

  ~OpenCoilSolver() override;

  void Init(hephaestus::GridFunctions & gridfunctions,
            const hephaestus::FESpaces & fespaces,
            hephaestus::BCMap & bc_map,
            hephaestus::Coefficients & coefficients) override;
  void Apply(mfem::ParLinearForm * lf) override;
  void SubtractSource(mfem::ParGridFunction * gf) override;

  // Initialises the child submesh.
  void initChildMesh();

  // Creates the relevant FE Collections and Spaces for the child submesh.
  void makeFESpaces();

  // Creates the relevant GridFunctions for the child submesh.
  void makeGridFunctions();

  // Sets up the boundary conditions to be used in the ScalarPotentialSource.
  // calculation.
  void setBCs();

  // Solves for the divergence-free Hodge dual of the electric current based on
  // Dirichlet BCs.
  void SPSCurrent();

  // Creates a mass matrix with basis functions that will be used in the Apply()
  // method
  void buildM1();

  // Sets the boundary attribute for the face to be used as reference in flux
  // calculation
  void setRefFace(const int face);

private:
  // Parameters
  int order_h1_;
  int order_hcurl_;
  int ref_face_;
  bool grad_phi_transfer_;

  std::pair<int, int> elec_attrs_;
  mfem::Array<int> coil_domains_;
  mfem::Array<int> coil_markers_;
  hephaestus::InputParameters solver_options_;

  mfem::Coefficient * sigma_{nullptr};
  mfem::Coefficient * Itotal_{nullptr};

  bool owns_sigma_{false};
  bool owns_Itotal_{false};

  // Names
  std::string grad_phi_name_;
  std::string V_gf_name_;
  std::string I_coef_name_;
  std::string cond_coef_name_;

  // Parent mesh, FE space, and current
  mfem::ParMesh * mesh_parent_{nullptr};

  mfem::ParGridFunction * grad_phi_parent_{nullptr};
  std::unique_ptr<mfem::ParGridFunction> grad_phi_t_parent_{nullptr};

  mfem::ParGridFunction * V_parent_{nullptr};
  std::unique_ptr<mfem::ParGridFunction> Vt_parent_{nullptr};

  // Child mesh and FE spaces
  std::unique_ptr<mfem::ParSubMesh> mesh_{nullptr};

  std::unique_ptr<mfem::ParFiniteElementSpace> H1FESpace_{nullptr};
  std::unique_ptr<mfem::H1_FECollection> H1FESpace_fec_{nullptr};

  std::unique_ptr<mfem::ParFiniteElementSpace> HCurlFESpace_{nullptr};
  std::unique_ptr<mfem::ND_FECollection> HCurlFESpace_fec{nullptr};

  // Child GridFunctions
  std::unique_ptr<mfem::ParGridFunction> grad_phi_{nullptr};
  std::unique_ptr<mfem::ParGridFunction> V_{nullptr};

  // Child boundary condition objects
  mfem::FunctionCoefficient high_src_;
  mfem::FunctionCoefficient low_src_;

  mfem::Array<int> high_terminal_;
  mfem::Array<int> low_terminal_;

  // Mass Matrix
  std::unique_ptr<mfem::ParBilinearForm> m1_{nullptr};

  // BC Map
  hephaestus::BCMap bc_maps;

  // Final LinearForm
  std::unique_ptr<mfem::ParLinearForm> final_lf_{nullptr};
};

} // namespace hephaestus