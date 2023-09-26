#pragma once
#include "source_base.hpp"

namespace hephaestus {

class ClosedCoilSolver : public hephaestus::Source {

public:

    ClosedCoilSolver(const hephaestus::InputParameters &params, 
                     const std::vector<hephaestus::Subdomain> &coil_dom,
                     const double J_total,
                     const int electrode_face);

    ~ClosedCoilSolver();

    void Init(hephaestus::GridFunctions &gridfunctions,
              const hephaestus::FESpaces &fespaces, 
              hephaestus::BCMap &bc_map,
              hephaestus::Coefficients &coefficients) override;
    void Apply(mfem::ParLinearForm *lf) override;
    void SubtractSource(mfem::ParGridFunction *gf) override;

    // Finds the electrode face and applies a single domain attribute to a 1-element layer adjacent to it.
    // Applies a different domain attribute to all other elements.
    // Also applies different boundary attributes on the two opposing faces of the layer, to act as Dirichlet BCs.
    void makeWedge();

    // Splits a full parallel mesh into two parallel submeshes by extracting a 1-element wide layer adjacent to the
    // electrode face.
    void initSubMeshes();

    // Checks if two faces share an edge.
    bool isAdjacent(const int f1, const int f2, const mfem::ParMesh* _par_mesh);

    void SubdomainVecToArray(const std::vector<hephaestus::Subdomain> &sd, mfem::Array<int> &arr);

    // Extracting a submesh often erases boundary attribute information. This method
    // ensures that boundary attributes are carried over from parent to child mesh
    void inheritBdrAttributes(const mfem::ParMesh* parent_mesh, mfem::ParSubMesh* child_mesh);

    private:

    // Parameters
    double Jtotal_;
    int elec_;
    std::vector<hephaestus::Subdomain> coil_domains_;

    // Grandparent mesh and FE spaces
    std::string hcurl_fespace_name;

    mfem::ParMesh* mesh_grandparent_;
    mfem::ParFiniteElementSpace* HCurlFESpace_grandparent_;

    // Parent mesh and FE spaces
    mfem::ParSubMesh* mesh_parent_;
    mfem::ParFiniteElementSpace* HCurlFESpace_parent_;

    // Children mesh and FE spaces
    std::vector<mfem::ParMesh*> mesh_;
    std::vector<mfem::ParFiniteElementSpace*> H1FESpace_;
    std::vector<mfem::ParFiniteElementSpace*> HCurlFESpace_;

};

} // namespace hephaestus