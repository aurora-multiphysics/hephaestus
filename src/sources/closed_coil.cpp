#include "closed_coil.hpp"

namespace hephaestus {

ClosedCoilSolver::ClosedCoilSolver(const hephaestus::InputParameters &params, 
                                   const std::vector<hephaestus::Subdomain> &coil_dom,
                                   const double Jtotal,
                                   const int electrode_face) 
                                   : hcurl_fespace_name(params.GetParam<std::string>("HCurlFESpaceName")),
                                     coil_domains_(coil_dom),
                                     Jtotal_(Jtotal),
                                     elec_(electrode_face) {

    std::cout << "Constructor" << std::endl;

}
// REMEMBER TO SET ALL POINTERS TO NULL

ClosedCoilSolver::~ClosedCoilSolver() {}

void ClosedCoilSolver::Init(hephaestus::GridFunctions &gridfunctions,
                           const hephaestus::FESpaces &fespaces, 
                           hephaestus::BCMap &bc_map,
                           hephaestus::Coefficients &coefficients) {

    std::cout << "Init" << std::endl;

    HCurlFESpace_grandparent_ = fespaces.Get(hcurl_fespace_name);
    if (HCurlFESpace_grandparent_ == NULL) {
        const std::string error_message = hcurl_fespace_name +
                                        " not found in fespaces when "
                                        "creating ClosedCoilSolver\n";
        mfem::mfem_error(error_message.c_str());
    }

    mesh_grandparent_ = HCurlFESpace_grandparent_->GetParMesh();

    std::cout << "Got mesh! Number of elements = " << mesh_grandparent_->GetNE() << std::endl;

    std::cout << "Preparing parent mesh..." << std::endl;
    mfem::Array<int> doms_array;
    SubdomainVecToArray(coil_domains_, doms_array);

    mesh_parent_ = new mfem::ParSubMesh(mfem::ParSubMesh::CreateFromDomain(*mesh_grandparent_,doms_array));

    inheritBdrAttributes(mesh_grandparent_, mesh_parent_);

    mesh_parent_->Save("test_parent.mesh");
    mesh_grandparent_->Save("test_grandparent.mesh");

    makeWedge();


   ////////////////// 


}

void ClosedCoilSolver::Apply(mfem::ParLinearForm *lf) {}

void ClosedCoilSolver::SubtractSource(mfem::ParGridFunction *gf) {}

void ClosedCoilSolver::SubdomainVecToArray(const std::vector<hephaestus::Subdomain> &sd, mfem::Array<int> &arr){
    
    arr.DeleteAll();
    for (auto s:sd) arr.Append(s.id);
}

void ClosedCoilSolver::inheritBdrAttributes(const mfem::ParMesh* parent_mesh, mfem::ParSubMesh* child_mesh){

    int face, ori, att;
    auto map = child_mesh->GetParentToSubMeshFaceIDMap();

    for (int bdr=0; bdr<parent_mesh->GetNBE(); ++bdr){

        parent_mesh->GetBdrElementFace(bdr, &face, &ori);
        if (map[face] != -1){
            att = parent_mesh->GetBdrAttribute(bdr);
            auto* new_elem = child_mesh->GetFace(map[face])->Duplicate(child_mesh);
            new_elem->SetAttribute(att);
            child_mesh->AddBdrElement(new_elem);
        }

    }

    child_mesh->FinalizeTopology();
    child_mesh->Finalize();
    child_mesh->SetAttributes();
}

void ClosedCoilSolver::makeWedge() {

    std::vector<int> bdr_els;
    int new_domain_attr = 2;
    std::cout << "Boundary attribute is " << elec_ << std::endl;

    // First we need to find the electrode boundary
    for (int i = 0; i<mesh_parent_->GetNBE(); ++i){
        if (mesh_parent_->GetBdrAttribute(i) == elec_){
            bdr_els.push_back(i);
            mesh_parent_->SetBdrAttribute(i, elec_);
        }
    }

    std::cout << "boundary elements = " << bdr_els.size() << std::endl;

    // Now we find all MPI ranks that are adjacent to some of the electrode and choose only the highest rank to contain the element-wide
    // transition region. This is explicitly for the case where the rank distribution cuts right through the electrode
    int bdr_thisrank = bdr_els.size();
    int *bdr_allranks = (int *)malloc(sizeof(int) * mfem::Mpi::WorldSize());
    MPI_Allgather(&bdr_thisrank, 1, MPI_INT, bdr_allranks, 1, MPI_INT, MPI_COMM_WORLD);

    int e_rank = 0;
    bool verify = false;
    for (int r=0; r<mfem::Mpi::WorldSize(); ++r) {
        if (bdr_allranks[r]) {
            e_rank = r;
            verify = true;
        }
    }

    if (!verify) MFEM_ABORT("Could not find electrode in any of the ranks!");
    free(bdr_allranks);

    // Now we make the wedge
    for (int e=0; e<mesh_parent_->GetNE(); ++e) mesh_parent_->SetAttribute(e, new_domain_attr-1);

    if (mfem::Mpi::WorldRank() == e_rank){

        int e1, e2;
        mfem::Array<int> e1_faces, ori;

        for (auto fc:bdr_els) {

            mesh_parent_->GetFaceElements(mesh_parent_->GetBdrFace(fc), &e1, &e2);
            mesh_parent_->SetAttribute(e1,new_domain_attr);

            // We need to find a face one element over and set its boundary attribute
            mesh_parent_->GetElementFaces(e1, e1_faces, ori);
            for (auto e1_fc:e1_faces) {
                if (e1_fc != mesh_parent_->GetBdrFace(fc) && !isAdjacent(e1_fc, mesh_parent_->GetBdrFace(fc), mesh_parent_)) {
                    auto* new_elem = mesh_parent_->GetFace(e1_fc)->Duplicate(mesh_parent_);
                    new_elem->SetAttribute(elec_+10);
                    mesh_parent_->AddBdrElement(new_elem);
                }
            }
        }
    }

    mesh_parent_->FinalizeTopology();
    mesh_parent_->Finalize();
    mesh_parent_->SetAttributes();

}

void ClosedCoilSolver::initSubMeshes() {

    
}

bool ClosedCoilSolver::isAdjacent(const int f1, const int f2, const mfem::ParMesh* _par_mesh) {

    mfem::Array<int> edges_f1, edges_f2, orientations;
    _par_mesh->GetFaceEdges(f1, edges_f1, orientations);
    _par_mesh->GetFaceEdges(f2, edges_f2, orientations);

    for (auto e1:edges_f1) {
        for (auto e2:edges_f2){
            if (e1 == e2) return true;    
        }
    }

    return false;

}

}; // namespace hephaestus