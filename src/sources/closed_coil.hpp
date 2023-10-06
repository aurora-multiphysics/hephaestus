#pragma once
#include "source_base.hpp"
#include "scalar_potential_source.hpp"

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

    // Checks if two faces share an edge.
    bool isAdjacent(const int f1, const int f2, const mfem::ParMesh* _par_mesh);

    // Takes in either a Subdomain or a vector of subdomains, and aranges them in an array for
    // creating submeshes
    void SubdomainToArray(const std::vector<hephaestus::Subdomain> &sd, mfem::Array<int> &arr);
    void SubdomainToArray(const hephaestus::Subdomain &sd, mfem::Array<int> &arr);

    // Takes in a boundary attribute and puts it into the kind of marker array used by MFEM
    void toMarkerArray(const int bdr_att, mfem::Array <int> &marker_array);

    // Takes in two marker arrays and returns the union of those
    void markerUnion(const mfem::Array<int> &m1, const mfem::Array<int> &m2, mfem::Array<int> &u);

    // Extracting a submesh often erases boundary attribute information. This method
    // ensures that boundary attributes are carried over from parent to child mesh
    void inheritBdrAttributes(const mfem::ParMesh* parent_mesh, mfem::ParSubMesh* child_mesh);

    // Solves for the divergence-free current based on Dirichlet BCs
    void solveLaplace();
    void SPSCurrent();

    // Creates the grad operator to calculate the current from the potential
    void buildGrad();

    // Creates the relevant Bilinear operators on the child submeshes
    void buildLaplace();

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

    // Sets the markers and DoFs for the Dirichlet BCs in the children submeshes 
    void setBCs();

    // Applies the gradient operator to the children submesh potentials, and performs a helmholz projection
    // to obtain divergence-free currents. Finally, normalises the current such that the total current
    // going through the electrode face is 1.
    void calcCurrent();

    // Applies a Helmholtz projection to a current to remove any divergences
    void removeJDivergence(mfem::ParGridFunction* Jraw, mfem::ParGridFunction* J);

    // Checks whether a given element is within a certain set of domains
    bool isInDomain(const int el, const std::vector<hephaestus::Subdomain> &dom, const mfem::ParMesh* mesh);
    bool isInDomain(const int el, const hephaestus::Subdomain &sd, const mfem::ParMesh* mesh);

    private:

    // Parameters
    double Jtotal_;
    int elec_;
    int order_;
    std::pair<int,int> elec_attrs_;
    std::vector<hephaestus::Subdomain> coil_domains_;
    int new_domain_attr_;
    mfem::ConstantCoefficient* coef1_;
    mfem::ConstantCoefficient* coef0_;

    // FE space, mesh, and J GridFunction
    std::string hcurl_fespace_name;
    std::string J_gf_name;

    // Parent mesh, FE space, and current
    mfem::ParMesh* mesh_parent_;
    mfem::ParFiniteElementSpace* HCurlFESpace_parent_;
    mfem::ParGridFunction* J_parent_;

    // Children mesh and FE spaces
    std::vector<mfem::ParSubMesh*> mesh_;
    std::vector<mfem::H1_FECollection*> H1_Collection_;
    std::vector<mfem::ND_FECollection*> HCurl_Collection_;
    std::vector<mfem::ParFiniteElementSpace*> H1FESpace_;
    std::vector<mfem::ParFiniteElementSpace*> HCurlFESpace_;

    // Children GridFunctions
    std::vector<mfem::ParGridFunction*> Jr_;
    std::vector<mfem::ParGridFunction*> J_;
    std::vector<mfem::ParGridFunction*> V_;

    // Children IR Order and geometry
    std::vector<int> irOrder_;
    std::vector<int> geom_;
    std::vector<mfem::IntegrationRule> intrule_;

    // Children boundaries and TDOFs
    std::vector<mfem::Array<int>> ess_bdr_;
    std::vector<mfem::Array<int>> ess_bdr_tdofs_;

    // Children Ops
    std::vector<mfem::ParDiscreteLinearOperator*> grad_;
    std::vector<mfem::ParBilinearForm*> laplace_;
    std::vector<mfem::ParLinearForm*> rhod_;

    // Children PCG objects
    std::vector<mfem::HypreParMatrix*> laplace_hypre_;
    std::vector<mfem::HypreParVector*> V_hypre_;
    std::vector<mfem::HypreParVector*> rhs_hypre_;
    std::vector<mfem::HypreBoomerAMG*> amg_;
    std::vector<mfem::HyprePCG*>       pcg_;

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