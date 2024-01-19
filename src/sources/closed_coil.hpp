#pragma once
#include "open_coil.hpp"
#include "scalar_potential_source.hpp"
#include "source_base.hpp"

namespace hephaestus
{

// Calculates the flux of a vector field v_field through the face with boundary
// attribute face_attr.
double calcFluxCC(mfem::GridFunction * v_field, int face_attr);

class ClosedCoilSolver : public hephaestus::Source
{

public:
  ClosedCoilSolver(const hephaestus::InputParameters & params,
                   const mfem::Array<int> & coil_dom,
                   const int electrode_face);

  // Override virtual Source destructor to avoid leaks.
  ~ClosedCoilSolver() override;

  void Init(hephaestus::GridFunctions & gridfunctions,
            const hephaestus::FESpaces & fespaces,
            hephaestus::BCMap & bc_map,
            hephaestus::Coefficients & coefficients) override;
  void Apply(mfem::ParLinearForm * lf) override;
  void SubtractSource(mfem::ParGridFunction * gf) override;

  // Finds the electrode face and applies a single domain attribute to a
  // 1-element layer adjacent to it. Applies a different domain attribute to
  // other elements in the coil. Also applies different boundary attributes on
  // the two opposing faces of the layer, to act as Dirichlet BCs.
  void makeWedge();

  // Extracts the coil submesh and prepares the gridfunctions and FE spaces
  // for being passed to the OpenCoilSolver in the transition region
  void prepareCoilSubmesh();

  // Applies the OpenCoilSolver to the transition region
  void solveTransition();

  // Solves for the current in the coil region
  void solveCoil();

  // Resets the domain attributes on the parent mesh to what they were initially
  void restoreAttributes();

  // Finds the coordinates for the "centre of mass" of the vertices of an
  // element.
  mfem::Vector elementCentre(int el, mfem::ParMesh * pm);

  // Checks whether a given element is within a certain domain or vector of
  // domains.
  bool isInDomain(const int el, const mfem::Array<int> & dom, const mfem::ParMesh * mesh);
  bool isInDomain(const int el, const int & sd, const mfem::ParMesh * mesh);

private:
  // Parameters
  int order_hcurl_;
  int order_h1_;
  int new_domain_attr_;
  std::pair<int, int> elec_attrs_;
  mfem::Array<int> coil_domains_;
  mfem::Array<int> coil_markers_;
  mfem::Array<int> transition_domain_;
  mfem::Array<int> transition_markers_;
  mfem::Coefficient * sigma_{nullptr};
  mfem::Coefficient * Itotal_{nullptr};
  std::vector<int> old_dom_attrs;
  hephaestus::InputParameters solver_options_;

  bool owns_sigma{false};
  bool owns_Itotal_{false};
  bool owns_grad_phi_parent_{false};
  bool owns_H1FESpace_parent_{false};

  // Here, we are solving for -(σ∇Va,∇ψ) = (σ∇Vt,∇ψ), where ∇Vt is grad_phi_t (within its relevant
  // mesh), ∇Va is grad_phi_aux, and their sum ∇Vt+∇Va is the full grad_phi, which is, up to an
  // overall sign, equal to the electric field in the electrostatic case. Seting grad_phi_transfer_
  // to true will negatively affect performance, but the final grad_phi function will be transferred
  // to a GridFunction for viewing purposes. Only set to true if you wish to visualise the final
  // grad_phi.
  bool grad_phi_transfer_{false};

  // Names
  std::string hcurl_fespace_name_;
  std::string cond_coef_name_;
  std::string h1_fespace_name_;
  std::string grad_phi_name_;
  std::string I_coef_name_;

  // Parent mesh, FE space, and current
  mfem::ParMesh * mesh_parent_{nullptr};
  mfem::ParGridFunction * grad_phi_parent_{nullptr};
  mfem::ParFiniteElementSpace * HCurlFESpace_parent_{nullptr};
  mfem::ParFiniteElementSpace * H1FESpace_parent_{nullptr};

  std::unique_ptr<mfem::H1_FECollection> H1FESpace_parent_fec_{nullptr};

  // In case J transfer is true
  std::unique_ptr<mfem::ParGridFunction> grad_phi_t_parent_{nullptr};

  // Coil mesh, FE Space, and current
  std::unique_ptr<mfem::ParSubMesh> mesh_coil_{nullptr};
  std::unique_ptr<mfem::ParSubMesh> mesh_t_{nullptr};
  std::unique_ptr<mfem::ParFiniteElementSpace> H1FESpace_coil_{nullptr};
  std::unique_ptr<mfem::ParGridFunction> grad_phi_aux_coil_{nullptr};
  std::unique_ptr<mfem::ParGridFunction> V_coil_{nullptr};

  std::unique_ptr<mfem::ND_FECollection> grad_phi_aux_coil_fec_{nullptr};
  std::unique_ptr<mfem::ParFiniteElementSpace> grad_phi_aux_coil_fes{nullptr};

  // Final LinearForm
  std::unique_ptr<mfem::ParLinearForm> final_lf_{nullptr};
};

class Plane3D
{

public:
  Plane3D();

  // Constructs a mathematical 3D plane from a mesh face
  void make3DPlane(const mfem::ParMesh * pm, const int face);

  // Calculates on which side of the infinite 3D plane a point is.
  // Returns 1, -1, or 0, the latter meaning the point is on the plane
  int side(const mfem::Vector v);

private:
  std::unique_ptr<mfem::Vector> u{nullptr};
  double d;
};

} // namespace hephaestus