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
                   mfem::Array<int> coil_dom,
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
  void MakeWedge();

  // Extracts the coil submesh and prepares the gridfunctions and FE spaces
  // for being passed to the OpenCoilSolver in the transition region
  void PrepareCoilSubmesh();

  // Applies the OpenCoilSolver to the transition region
  void SolveTransition();

  // Solves for the current in the coil region
  void SolveCoil();

  // Resets the domain attributes on the parent mesh to what they were initially
  void RestoreAttributes();

  // Finds the coordinates for the "centre of mass" of the vertices of an
  // element.
  mfem::Vector ElementCentre(int el, mfem::ParMesh * pm);

  // Checks whether a given element is within a certain domain or vector of
  // domains.
  bool IsInDomain(const int el, const mfem::Array<int> & dom, const mfem::ParMesh * mesh);
  bool IsInDomain(const int el, const int & sd, const mfem::ParMesh * mesh);

private:
  // Parameters
  int _order_hcurl;
  int _order_h1;
  int _new_domain_attr;
  std::pair<int, int> _elec_attrs;
  mfem::Array<int> _coil_domains;
  mfem::Array<int> _coil_markers;
  mfem::Array<int> _transition_domain;
  mfem::Array<int> _transition_markers;
  mfem::Coefficient * _sigma{nullptr};
  mfem::Coefficient * _itotal{nullptr};
  std::vector<int> _old_dom_attrs;
  hephaestus::InputParameters _solver_options;

  bool _owns_sigma{false};
  bool _owns_itotal{false};
  bool _owns_grad_phi_parent{false};
  bool _owns_h1_fe_space_parent{false};

  // Here, we are solving for -(σ∇Va,∇ψ) = (σ∇Vt,∇ψ), where ∇Vt is grad_phi_t (within its relevant
  // mesh), ∇Va is grad_phi_aux, and their sum ∇Vt+∇Va is the full grad_phi, which is, up to an
  // overall sign, equal to the electric field in the electrostatic case. Seting grad_phi_transfer_
  // to true will negatively affect performance, but the final grad_phi function will be transferred
  // to a GridFunction for viewing purposes. Only set to true if you wish to visualise the final
  // grad_phi.
  bool _grad_phi_transfer{false};

  // Names
  std::string _hcurl_fespace_name;
  std::string _cond_coef_name;
  std::string _h1_fespace_name;
  std::string _grad_phi_name;
  std::string _i_coef_name;

  // Parent mesh, FE space, and current
  mfem::ParMesh * _mesh_parent{nullptr};
  mfem::ParGridFunction * _grad_phi_parent{nullptr};
  mfem::ParFiniteElementSpace * _h_curl_fe_space_parent{nullptr};
  mfem::ParFiniteElementSpace * _h1_fe_space_parent{nullptr};

  std::unique_ptr<mfem::H1_FECollection> _h1_fe_space_parent_fec{nullptr};

  // In case J transfer is true
  std::unique_ptr<mfem::ParGridFunction> _grad_phi_t_parent{nullptr};

  // Coil mesh, FE Space, and current
  std::unique_ptr<mfem::ParSubMesh> _mesh_coil{nullptr};
  std::unique_ptr<mfem::ParSubMesh> _mesh_t{nullptr};
  std::unique_ptr<mfem::ParFiniteElementSpace> _h1_fe_space_coil{nullptr};
  std::unique_ptr<mfem::ParGridFunction> _grad_phi_aux_coil{nullptr};
  std::unique_ptr<mfem::ParGridFunction> _v_coil{nullptr};

  std::unique_ptr<mfem::ND_FECollection> _grad_phi_aux_coil_fec{nullptr};
  std::unique_ptr<mfem::ParFiniteElementSpace> _grad_phi_aux_coil_fes{nullptr};

  // Final LinearForm
  std::unique_ptr<mfem::ParLinearForm> _final_lf{nullptr};
};

class Plane3D
{

public:
  Plane3D();

  // Constructs a mathematical 3D plane from a mesh face
  void Make3DPlane(const mfem::ParMesh * pm, const int face);

  // Calculates on which side of the infinite 3D plane a point is.
  // Returns 1, -1, or 0, the latter meaning the point is on the plane
  int Side(const mfem::Vector v);

private:
  std::unique_ptr<mfem::Vector> _u{nullptr};
  double _d{0};
};

} // namespace hephaestus