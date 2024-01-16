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

template <typename T>
void ifDelete(T * ptr);

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
                 const mfem::Array<int> & coil_dom,
                 const std::pair<int, int> electrodes);

  ~OpenCoilSolver();

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
  mfem::Coefficient *sigma_;
  mfem::Coefficient *Itotal_;
  hephaestus::InputParameters solver_options_;

  // Names
  std::string grad_phi_name_;
  std::string V_gf_name_;
  std::string I_coef_name_;
  std::string cond_coef_name_;

  // Parent mesh, FE space, and current
  mfem::ParMesh *mesh_parent_;
  mfem::ParGridFunction *grad_p_parent_;
  mfem::ParGridFunction *grad_p_t_parent_;
  mfem::ParGridFunction *V_parent_;
  mfem::ParGridFunction *Vt_parent_;

  // Child mesh and FE spaces
  mfem::ParSubMesh * mesh_;
  mfem::ParFiniteElementSpace * H1FESpace_;
  mfem::ParFiniteElementSpace * HCurlFESpace_;

  // Child GridFunctions
  mfem::ParGridFunction *grad_p_;
  mfem::ParGridFunction *V_;

  // Child boundary condition objects
  mfem::FunctionCoefficient high_src_;
  mfem::FunctionCoefficient low_src_;
  mfem::Array<int> high_terminal_;
  mfem::Array<int> low_terminal_;

  // Mass Matrix
  mfem::ParBilinearForm * m1_;

  // BC Map
  hephaestus::BCMap bc_maps;

  // Final LinearForm
  mfem::ParLinearForm * final_lf_;
};

} // namespace hephaestus