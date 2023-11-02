#pragma once
#include "helmholtz_projector.hpp"
#include "scalar_potential_source.hpp"
#include "source_base.hpp"

namespace hephaestus {

double calcFlux(mfem::GridFunction *v_field, int face_attr);

template <typename T> void ifDelete(T *ptr);

void inheritBdrAttributes(const mfem::ParMesh *parent_mesh,
                          mfem::ParSubMesh *child_mesh);

class OpenCoilSolver : public hephaestus::Source {

public:
  OpenCoilSolver(const hephaestus::InputParameters &params,
                 const mfem::Array<int> &coil_dom,
                 const std::pair<int, int> electrodes, const int order);

  ~OpenCoilSolver();

  void Init(hephaestus::GridFunctions &gridfunctions,
            const hephaestus::FESpaces &fespaces, hephaestus::BCMap &bc_map,
            hephaestus::Coefficients &coefficients) override;
  void Apply(mfem::ParLinearForm *lf) override;
  void SubtractSource(mfem::ParGridFunction *gf) override;

  // Initialises the child submesh.
  void initChildMesh();

  // Creates the relevant FE Collections and Spaces for the child submesh.
  void makeFESpaces();

  // Creates the relevant GridFunctions for the child submesh.
  void makeGridFunctions();

  // Sets up the boundary conditions to be used in the ScalarPotentialSource.
  // calculation.
  void setBCs();

  // Solves for the divergence-free current based on Dirichlet BCs.
  void SPSCurrent();

  // Applies the HelmholtzProjector onto the J GridFunction to clean it of any
  // divergences
  void cleanDivergence(hephaestus::GridFunctions *gridfunctions,
                       std::string J_name, std::string V_name,
                       hephaestus::BCMap *bc_map);

  // Sets the boundary attribute for the face to be used as reference in flux
  // calculation
  void setRefFace(const int face);

private:
  // Parameters
  int order_;
  int ref_face_;
  std::pair<int, int> elec_attrs_;
  mfem::Array<int> coil_domains_;
  mfem::ConstantCoefficient *coef1_;
  mfem::ConstantCoefficient *coef0_;
  mfem::Coefficient *Itotal_;

  // Names
  std::string J_gf_name_;
  std::string V_gf_name_;
  std::string I_coef_name_;

  // Parent mesh, FE space, and current
  mfem::ParMesh *mesh_parent_;
  mfem::ParGridFunction *J_parent_;
  mfem::ParGridFunction *V_parent_;

  // Child mesh and FE spaces
  mfem::ParSubMesh *mesh_;
  mfem::H1_FECollection *H1_Collection_;
  mfem::ND_FECollection *HCurl_Collection_;
  mfem::ParFiniteElementSpace *H1FESpace_;
  mfem::ParFiniteElementSpace *HCurlFESpace_;

  // Child GridFunctions
  mfem::ParGridFunction *J_;
  mfem::ParGridFunction *V_;

  // Child boundary condition objects
  mfem::FunctionCoefficient *high_src_;
  mfem::FunctionCoefficient *low_src_;
  mfem::Array<int> high_terminal_;
  mfem::Array<int> low_terminal_;

  // Child ScalarPotentialSource objects
  hephaestus::ScalarPotentialSource *sps_;
  hephaestus::InputParameters *sps_params_;
  hephaestus::InputParameters *current_solver_options_;

  hephaestus::GridFunctions *gridfunctions_;
  hephaestus::FESpaces *fespaces_;
  hephaestus::BCMap *bc_maps_;
  hephaestus::Coefficients *coefs_;
  hephaestus::FunctionDirichletBC *high_DBC_;
  hephaestus::FunctionDirichletBC *low_DBC_;
};

} // namespace hephaestus