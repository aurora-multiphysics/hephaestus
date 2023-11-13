#pragma once
#include "open_coil.hpp"
#include "scalar_potential_source.hpp"
#include "source_base.hpp"

namespace hephaestus {

// Calculates the flux of a vector field v_field through the face with boundary
// attribute face_attr.
double calcFluxCC(mfem::GridFunction *v_field, int face_attr);

class ClosedCoilSolver : public hephaestus::Source {

public:
  ClosedCoilSolver(const hephaestus::InputParameters &params,
                   const mfem::Array<int> &coil_dom, const int electrode_face);

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

  // Extracts the coil submesh and prepares the gridfunctions and FE spaces
  // for being passed to the OpenCoilSolver in the transition region
  void prepareCoilSubmesh();

  // Applies the OpenCoilSolver to the transition region
  void solveTransition();

  // Solves for the current in the coil region
  void solveCoil();

  // Creates a mass matrix with basis functions that will be used in the Apply()
  // method
  void buildM1();

  // Sets the current such that its flux across the electrode face is 1
  void normaliseCurrent();

  // Resets the domain attributes on the parent mesh to what they were initially
  void restoreAttributes();

  // Finds the coordinates for the "centre of mass" of the vertices of an
  // element.
  mfem::Vector elementCentre(int el, mfem::ParMesh *pm);

  // Checks whether a given element is within a certain domain or vector of
  // domains.
  bool isInDomain(const int el, const mfem::Array<int> &dom,
                  const mfem::ParMesh *mesh);
  bool isInDomain(const int el, const int &sd, const mfem::ParMesh *mesh);

private:
  // Parameters
  int order_hcurl_;
  int order_h1_;
  int new_domain_attr_;
  std::pair<int, int> elec_attrs_;
  mfem::Array<int> coil_domains_;
  mfem::Array<int> transition_domain_;
  mfem::Coefficient *Itotal_;
  std::vector<int> old_dom_attrs;

  // Names
  std::string hcurl_fespace_name_;
  std::string J_gf_name_;
  std::string I_coef_name_;

  // Parent mesh, FE space, and current
  mfem::ParMesh *mesh_parent_;
  mfem::ParGridFunction *J_parent_;
  mfem::ParFiniteElementSpace *HCurlFESpace_parent_;

  // Coil mesh, FE Space, and current
  mfem::ParSubMesh *mesh_coil_;
  mfem::ParFiniteElementSpace *HCurlFESpace_coil_;
  mfem::ParFiniteElementSpace *H1FESpace_coil_;
  mfem::ParGridFunction *J_coil_;
  mfem::ParGridFunction *Jt_coil_;

  // Mass Matrix
  mfem::ParBilinearForm *m1_;
};

class Plane3D {

public:
  Plane3D();
  ~Plane3D();

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