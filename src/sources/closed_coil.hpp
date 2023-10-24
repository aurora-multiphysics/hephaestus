#pragma once
#include "helmholtz_projector.hpp"
#include "div_free_source.hpp"
#include "scalar_potential_source.hpp"
#include "source_base.hpp"

namespace hephaestus {

// Calculates the flux of a vector field v_field through the face with boundary
// attribute face_attr.
double calcFlux(mfem::GridFunction *v_field, int face_attr);

class ClosedCoilSolver : public hephaestus::Source {

public:
  ClosedCoilSolver(const hephaestus::InputParameters &params,
                   const std::vector<hephaestus::Subdomain> &coil_dom,
                   const int electrode_face, const int order);

  ~ClosedCoilSolver();

  void Init(hephaestus::GridFunctions &gridfunctions,
            const hephaestus::FESpaces &fespaces, hephaestus::BCMap &bc_map,
            hephaestus::Coefficients &coefficients) override;
  void Apply(mfem::ParLinearForm *lf) override;
  void SubtractSource(mfem::ParGridFunction *gf) override;

  // Finds the electrode face and applies a single domain attribute to a
  // 1-element layer adjacent to it. Applies a different domain attribute to
  // other elements in the coil. Also applies different boundary attributes on
  // the two opposing faces of the layer, to act as Dirichlet BCs.
  void makeWedge();

  // Splits a full parallel mesh into two parallel submeshes by extracting a
  // 1-element wide layer adjacent to the electrode face as one submesh and the
  // rest of the coil as the second submesh.
  void initChildMeshes();

  // Takes in either a Subdomain or a vector of subdomains, and arranges them in
  // an array for creating submeshes.
  void SubdomainToArray(const std::vector<hephaestus::Subdomain> &sd,
                        mfem::Array<int> &arr);
  void SubdomainToArray(const hephaestus::Subdomain &sd, mfem::Array<int> &arr);

  // Extracting a submesh sometimes erases boundary attribute information. This
  // method ensures that boundary attributes are carried over from parent to
  // child mesh.
  void inheritBdrAttributes(const mfem::ParMesh *parent_mesh,
                            mfem::ParSubMesh *child_mesh);

  // Solves for the divergence-free current based on Dirichlet BCs.
  void SPSCurrent();

  // Finds the coordinates for the "centre of mass" of the vertices of an
  // element.
  mfem::Vector elementCentre(int el, mfem::ParMesh *pm);

  // Resizes all of the vectors to the number of children submeshes (2, in this
  // case).
  void resizeChildVectors();

  // Creates the relevant FE Collections and Spaces for the child submeshes.
  void makeFESpaces();

  // Creates the relevant Grid Functions for the children submeshes.
  void makeGridFunctions();

  // Checks whether a given element is within a certain domain or vector of
  // domains.
  bool isInDomain(const int el, const std::vector<hephaestus::Subdomain> &dom,
                  const mfem::ParMesh *mesh);
  bool isInDomain(const int el, const hephaestus::Subdomain &sd,
                  const mfem::ParMesh *mesh);

  // Sets up the boundary conditions to be used in the ScalarPotentialSource
  // calculation.
  void SetBCs();

  // Resets the domain attributes on the parent mesh to what they were initially
  void restoreAttributes();

  // Applies the HelmholtzProjector onto the J GridFunction to clean it of any divergences 
  void cleanDivergence(hephaestus::GridFunctions* gridfunctions, std::string J_name, std::string V_name, hephaestus::BCMap* bc_map);

private:
  // Parameters
  int order_;
  int new_domain_attr_;
  std::pair<int, int> elec_attrs_;
  std::vector<hephaestus::Subdomain> coil_domains_;
  mfem::ConstantCoefficient *coef1_;
  mfem::ConstantCoefficient *coef0_;
  mfem::Coefficient *Itotal_;
  std::vector<int> old_dom_attrs;

  // FE space, mesh, and J GridFunction
  std::string hcurl_fespace_name_;
  std::string J_gf_name_;
  std::string I_coef_name_;

  // Parent mesh, FE space, and current
  mfem::ParMesh *mesh_parent_;
  mfem::ParGridFunction *J_parent_;
  mfem::ParFiniteElementSpace *HCurlFESpace_parent_;

  // Children mesh and FE spaces
  std::vector<mfem::ParSubMesh *> mesh_;
  std::vector<mfem::H1_FECollection *> H1_Collection_;
  std::vector<mfem::ND_FECollection *> HCurl_Collection_;
  std::vector<mfem::ParFiniteElementSpace *> H1FESpace_;
  std::vector<mfem::ParFiniteElementSpace *> HCurlFESpace_;

  // Children GridFunctions
  std::vector<mfem::ParGridFunction *> J_;
  std::vector<mfem::ParGridFunction *> V_;

  // Children ScalarPotentialSource objects
  std::vector<hephaestus::InputParameters *> sps_params_;
  std::vector<hephaestus::InputParameters *> current_solver_options_;
  std::vector<hephaestus::ScalarPotentialSource *> sps_;
  std::vector<hephaestus::GridFunctions *> gridfunctions_;
  std::vector<hephaestus::FESpaces *> fespaces_;
  std::vector<hephaestus::BCMap *> bc_maps_;
  std::vector<hephaestus::Coefficients *> coefs_;
  std::vector<hephaestus::FunctionDirichletBC *> high_DBC_;
  std::vector<hephaestus::FunctionDirichletBC *> low_DBC_;

  // Children boundary condition objects
  mfem::FunctionCoefficient *high_src_;
  mfem::FunctionCoefficient *low_src_;
  mfem::VectorGridFunctionCoefficient *src_J_;

  mfem::Array<int> high_terminal_;
  mfem::Array<int> low_terminal_;
};

class Plane3D {

public:
  Plane3D();
  ~Plane3D(){};

  // Constructs a mathematical 3D plane from a mesh face
  void make3DPlane(const mfem::ParMesh *pm, const int face);

  // Calculates on which side of the infinite 3D plane a point is.
  // Returns 1, -1, or 0, the latter meaning the point is on the plane
  int side(const mfem::Vector v);

private:
  mfem::Vector *u;
  double d;
};

} // namespace hephaestus