#pragma once
#include "source_base.hpp"
#include "scalar_potential_source.hpp"
#include "div_free_source.hpp"

namespace hephaestus {

class ClosedCoilSolver : public hephaestus::Source {

public:

    ClosedCoilSolver(const hephaestus::InputParameters &params, 
                     const std::vector<hephaestus::Subdomain> &coil_dom,
                     const double J_total,
                     const int electrode_face,
                     const int order);

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
    void initChildMeshes();

    // Takes in either a Subdomain or a vector of subdomains, and aranges them in an array for
    // creating submeshes
    void SubdomainToArray(const std::vector<hephaestus::Subdomain> &sd, mfem::Array<int> &arr);
    void SubdomainToArray(const hephaestus::Subdomain &sd, mfem::Array<int> &arr);

    // Extracting a submesh often erases boundary attribute information. This method
    // ensures that boundary attributes are carried over from parent to child mesh
    void inheritBdrAttributes(const mfem::ParMesh* parent_mesh, mfem::ParSubMesh* child_mesh);

    // Solves for the divergence-free current based on Dirichlet BCs
    void SPSCurrent();

    // Calculates the flux of J through the face with attribute face_attr within the child submesh
    // with index idx
    double calcJFlux(int face_attr, int idx);

    // Finds the coordinates for the "centre of mass" of the vertices of an element
    mfem::Vector elementCentre(int el, mfem::ParMesh* pm);

    // Resizes all of the vectors to the number of children submeshes (2, in this case)
    void resizeChildVectors();

    // Creates the relevant FE Collections and Spaces for the children submeshes
    void makeFESpaces();

    // Creates the relevant Grid Functions for the children submeshes
    void makeGridFunctions();

    // Applies the gradient operator to the children submesh potentials, and performs a helmholz projection
    // to obtain divergence-free currents. Finally, normalises the current such that the total current
    // going through the electrode face is 1.
    void calcCurrent();

    // Checks whether a given element is within a certain set of domains
    bool isInDomain(const int el, const std::vector<hephaestus::Subdomain> &dom, const mfem::ParMesh* mesh);
    bool isInDomain(const int el, const hephaestus::Subdomain &sd, const mfem::ParMesh* mesh);

    // Initialises the pointers for SPS parameters and current solver options
    void initSPSOptions();

    // Sets up the boundary conditions to be used in the ScalarPotentialSource calculation
    void SetBCs();

    private:

    // Parameters
    double Jtotal_;
    int order_;
    int new_domain_attr_;
    std::pair<int,int> elec_attrs_;
    std::vector<hephaestus::Subdomain> coil_domains_;
    mfem::ConstantCoefficient* coef1_;
    mfem::ConstantCoefficient* coef0_;

    // FE space, mesh, and J GridFunction
    std::string hcurl_fespace_name;
    std::string J_gf_name;

    // Parent mesh, FE space, and current
    mfem::ParMesh* mesh_parent_;
    mfem::ParFiniteElementSpace* HCurlFESpace_parent_;
    mfem::ParGridFunction* J_parent_;
    mfem::ParGridFunction* Jr_parent_;

    // Children mesh and FE spaces
    std::vector<mfem::ParSubMesh*> mesh_;
    std::vector<mfem::H1_FECollection*> H1_Collection_;
    std::vector<mfem::ND_FECollection*> HCurl_Collection_;
    std::vector<mfem::ParFiniteElementSpace*> H1FESpace_;
    std::vector<mfem::ParFiniteElementSpace*> HCurlFESpace_;

    // Children GridFunctions
    std::vector<mfem::ParGridFunction*> J_;
    std::vector<mfem::ParGridFunction*> V_;

    // Children ScalarPotentialSource objects
    std::vector<hephaestus::InputParameters*> sps_params_;
    std::vector<hephaestus::InputParameters*> current_solver_options_;
    std::vector<hephaestus::ScalarPotentialSource*> sps_;
    std::vector<hephaestus::GridFunctions*> gridfunctions_;
    std::vector<hephaestus::FESpaces*> fespaces_;
    std::vector<hephaestus::BCMap*> bc_maps_;
    std::vector<hephaestus::Coefficients*> coefs_;
    std::vector<hephaestus::FunctionDirichletBC*> high_DBC_;
    std::vector<hephaestus::FunctionDirichletBC*> low_DBC_;

    mfem::FunctionCoefficient* high_src_;
    mfem::FunctionCoefficient* low_src_;
    mfem::Array<int> high_terminal_;
    mfem::Array<int> low_terminal_;


};


class Plane3D {

    public:

    
    Plane3D(); 
    ~Plane3D(){};

    // Constructs a mathematical 3D plane from a mesh face
    void make3DPlane(const mfem::ParMesh* pm, const int face);

    // Calculates on which side of the infinite 3D plane a point is.
    // Returns 1, -1, or 0, the latter meaning the point is on the plane
    int side(const mfem::Vector v);

    private:

    mfem::Vector* u;
    double d;

};

} // namespace hephaestus