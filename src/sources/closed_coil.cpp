#include "closed_coil.hpp"

namespace hephaestus {

ClosedCoilSolver::ClosedCoilSolver(const hephaestus::InputParameters &params, 
                                   const std::vector<hephaestus::Subdomain> &coil_dom,
                                   const double Jtotal,
                                   const int electrode_face,
                                   const int order) 
                                   : hcurl_fespace_name(params.GetParam<std::string>("HCurlFESpaceName")),
                                     J_gf_name(params.GetParam<std::string>("JGridFunctionName")),
                                     coil_domains_(coil_dom),
                                     Jtotal_(Jtotal),
                                     elec_(electrode_face),
                                     order_(order) {

    std::cout << "Constructor" << std::endl;

}
// REMEMBER TO SET ALL POINTERS TO NULL

ClosedCoilSolver::~ClosedCoilSolver() {}

void ClosedCoilSolver::Init(hephaestus::GridFunctions &gridfunctions,
                           const hephaestus::FESpaces &fespaces, 
                           hephaestus::BCMap &bc_map,
                           hephaestus::Coefficients &coefficients) {

    std::cout << "Init" << std::endl;
    coef1_ = new mfem::ConstantCoefficient(1.0);
    coef0_ = new mfem::ConstantCoefficient(0.0);

    // Retrieving the grandparent FE space and mesh
    HCurlFESpace_grandparent_ = fespaces.Get(hcurl_fespace_name);
    if (HCurlFESpace_grandparent_ == NULL) {
        const std::string error_message = hcurl_fespace_name +
                                        " not found in fespaces when "
                                        "creating ClosedCoilSolver\n";
        mfem::mfem_error(error_message.c_str());
    }

    J_grandparent_ = gridfunctions.Get(J_gf_name);
    if (J_grandparent_ == NULL) {
        const std::string error_message = J_gf_name +
                                        " not found in gridfunctions when "
                                        "creating ClosedCoilSolver\n";
        mfem::mfem_error(error_message.c_str());
    }


    mesh_grandparent_ = HCurlFESpace_grandparent_->GetParMesh();

    // Extracting parent submesh and setting FE spaces and grid functions
    mfem::Array<int> doms_array;
    SubdomainToArray(coil_domains_, doms_array);

    mesh_parent_ = new mfem::ParSubMesh(mfem::ParSubMesh::CreateFromDomain(*mesh_grandparent_,doms_array));
    inheritBdrAttributes(mesh_grandparent_, mesh_parent_);

    HCurl_Collection_parent_ = new mfem::ND_FECollection(order_, mesh_parent_->Dimension());
    HCurlFESpace_parent_ = new mfem::ParFiniteElementSpace(mesh_parent_, HCurl_Collection_parent_);
    J_ref_parent_ = new mfem::ParGridFunction(HCurlFESpace_parent_);
    *J_ref_parent_ = 0.0;

    // Now we extract the child submeshes
    makeWedge();
    initChildMeshes();
    solveCurrent();

    for (int i=0; i<2; ++i) mesh_[i]->Transfer(*J_[i], *J_ref_parent_); 

    //////////////////

    mesh_grandparent_->Save("test_grandparent.mesh");
    mesh_parent_->Save("test_parent.mesh");
    mesh_[0]->Save("test_child_1.mesh");
    mesh_[1]->Save("test_child_2.mesh");
 
}

void ClosedCoilSolver::Apply(mfem::ParLinearForm *lf) {

    mesh_parent_->Transfer(*J_ref_parent_, *J_grandparent_);
    *J_grandparent_ *= Jtotal_;
}

void ClosedCoilSolver::SubtractSource(mfem::ParGridFunction *gf) {}

void ClosedCoilSolver::SubdomainToArray(const std::vector<hephaestus::Subdomain> &sd, mfem::Array<int> &arr){
    
    arr.DeleteAll();
    for (auto s:sd) arr.Append(s.id);
}

void ClosedCoilSolver::SubdomainToArray(const hephaestus::Subdomain &sd, mfem::Array<int> &arr){
    
    arr.DeleteAll();
    arr.Append(sd.id);
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
    elec_attrs_.first = elec_;
    elec_attrs_.second = mesh_parent_->bdr_attributes.Max()+1;

    // First we need to find the electrode boundary
    for (int i = 0; i<mesh_parent_->GetNBE(); ++i){
        if (mesh_parent_->GetBdrAttribute(i) == elec_attrs_.first){
            bdr_els.push_back(i);
            //mesh_parent_->SetBdrAttribute(i, elec_attrs_.first);
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
                    new_elem->SetAttribute(elec_attrs_.second);
                    mesh_parent_->AddBdrElement(new_elem);
                }
            }
        }
    }

    mesh_parent_->FinalizeTopology();
    mesh_parent_->Finalize();
    mesh_parent_->SetAttributes();

}

void ClosedCoilSolver::initChildMeshes() {

    mfem::Array<int> doms_array;

    mesh_.resize(2);
    for (int i=0; i<2; ++i){
        hephaestus::Subdomain wedge("wedge", i+1);
        SubdomainToArray(wedge, doms_array);
        mesh_[i] = new mfem::ParSubMesh(mfem::ParSubMesh::CreateFromDomain(*mesh_parent_,doms_array));
        inheritBdrAttributes(mesh_parent_, mesh_[i]);
    }

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

void ClosedCoilSolver::toMarkerArray(const int bdr_att, mfem::Array <int> &marker_array){

    marker_array = 0;
    marker_array[bdr_att-1] = 1;    
}

void ClosedCoilSolver::markerUnion(const mfem::Array<int> &m1, const mfem::Array<int> &m2, mfem::Array<int> &u){

    assert(m1.Size() == m2.Size());
    int size = m1.Size();
    u.SetSize(size);
    for (int k=0; k<size; k++) {
        if ((m1[k] == 0) && (m2[k] == 0)) u[k] = 0;
        else u[k] = 1;
    }
}

void ClosedCoilSolver::solveCurrent(){

    // THIS IS MORE OF A TEMPORARY FUNCTION
    // IN THE FUTURE, IMPLEMENT THIS AS A ScalarPotentialSource
    H1_Collection_.resize(2);
    HCurl_Collection_.resize(2);
    H1FESpace_.resize(2);
    HCurlFESpace_.resize(2);

    irOrder_.resize(2);
    geom_.resize(2);
    intrule_.resize(2);

    J_.resize(2);
    Jr_.resize(2);
    V_.resize(2);

    ess_bdr_.resize(2);
    ess_bdr_tdofs_.resize(2);

    grad_.resize(2);
    laplace_.resize(2);
    rhod_.resize(2);

    laplace_hypre_.resize(2);
    V_hypre_.resize(2);
    rhs_hypre_.resize(2);
    amg_.resize(2);
    pcg_.resize(2);
    
    for (int i=0; i<2; ++i){

        std::cout << "Submesh " << i << ", FE Spaces" << std::endl;

        // FE spaces and grid functions
        H1_Collection_[i] = new mfem::H1_FECollection(order_, mesh_[i]->Dimension());
        HCurl_Collection_[i] = new mfem::ND_FECollection(order_, mesh_[i]->Dimension());
        H1FESpace_[i] = new mfem::ParFiniteElementSpace(mesh_[i], H1_Collection_[i]);
        HCurlFESpace_[i] = new mfem::ParFiniteElementSpace(mesh_[i], HCurl_Collection_[i]);

        std::cout << "Submesh " << i << ", GridFunctions" << std::endl;

        V_[i] = new mfem::ParGridFunction(H1FESpace_[i]);
        Jr_[i] = new mfem::ParGridFunction(HCurlFESpace_[i]);
        J_[i] = new mfem::ParGridFunction(HCurlFESpace_[i]);

        *V_[i] = 0.0;
        *Jr_[i] = 0.0;
        *J_[i] = 0.0;

        if (mesh_[i]->GetNE() != 0){
            irOrder_[i] = H1FESpace_[i]->GetElementTransformation(0)->OrderW() + 2*order_;
            geom_[i] = H1FESpace_[i]->GetFE(0)->GetGeomType();
            intrule_[i] = mfem::IntRules.Get(geom_[i],irOrder_[i]);
        }

        std::cout << "Submesh " << i << ", BCs" << std::endl;

        // Dirichlet BCs
        mfem::Array<int> bdr_array_V1(mesh_[i]->bdr_attributes.Max());
        mfem::Array<int> bdr_array_V0(mesh_[i]->bdr_attributes.Max());

        toMarkerArray(elec_attrs_.first, bdr_array_V0);
        toMarkerArray(elec_attrs_.second, bdr_array_V1);

        V_[i]->ProjectBdrCoefficient(*coef0_, bdr_array_V0);
        V_[i]->ProjectBdrCoefficient(*coef1_, bdr_array_V1);

        markerUnion(bdr_array_V0, bdr_array_V1, ess_bdr_[i]);
        H1FESpace_[i]->GetEssentialTrueDofs(ess_bdr_[i], ess_bdr_tdofs_[i]);
        
        std::cout << "Submesh " << i << ", Ops" << std::endl;

        // Setting up ops
        grad_[i] = new mfem::common::ParDiscreteGradOperator(H1FESpace_[i],HCurlFESpace_[i]);
        laplace_[i] = new mfem::ParBilinearForm(H1FESpace_[i]);
        rhod_[i] = new mfem::ParLinearForm(H1FESpace_[i]);

        *rhod_[i] = 0.0;

        std::cout << "Submesh " << i << ", Add integrators" << std::endl;


        laplace_[i]->AddDomainIntegrator(new mfem::DiffusionIntegrator(*coef1_));
        rhod_[i]->AddDomainIntegrator(new mfem::DomainLFIntegrator(*coef0_));

        std::cout << "Submesh " << i << ", Assemble" << std::endl;

        // Assembling everything
        rhod_[i]->Assemble();
        grad_[i]->Assemble();
        grad_[i]->Finalize();
        laplace_[i]->Assemble();
        laplace_[i]->Finalize();

        std::cout << "Submesh " << i << ", Solve" << std::endl;

        // At last, we solve the system using PCG
        laplace_hypre_[i] = new mfem::HypreParMatrix;
        V_hypre_[i] = new mfem::HypreParVector(H1FESpace_[i]);
        rhs_hypre_[i] = new mfem::HypreParVector(H1FESpace_[i]);
        laplace_[i]->FormLinearSystem(ess_bdr_tdofs_[i], *V_[i], *rhod_[i], 
                                        *laplace_hypre_[i], *V_hypre_[i], *rhs_hypre_[i]);
        amg_[i] = new mfem::HypreBoomerAMG(*laplace_hypre_[i]);
        pcg_[i] = new mfem::HyprePCG(*laplace_hypre_[i]);
        pcg_[i]->SetTol(1e-12);
        pcg_[i]->SetMaxIter(500);
        pcg_[i]->SetPrintLevel(2);
        pcg_[i]->SetPreconditioner(*amg_[i]);
        pcg_[i]->Mult(*rhs_hypre_[i], *V_hypre_[i]);
        laplace_[i]->RecoverFEMSolution(*V_hypre_[i],*rhod_[i],*V_[i]);     
        grad_[i]->Mult(*V_[i],*Jr_[i]);

        // THIS NEEDS TO BE MADE INTO A DIVERGENCELESS CURRENT
        *J_[i] = *Jr_[i];

        double flux = calcJFlux(elec_attrs_.first, i);

        *J_[i] /= abs(flux);
    }

    *J_[1] *= -1.0;

}

double ClosedCoilSolver::calcJFlux(int face_attr, int idx){

    double flux = 0.0;
    double area = 0.0;

    MFEM_ASSERT(HCurlFESpace_[idx].GetVDim() == 1, "");
    mfem::Vector local_dofs, normal_vec;
    mfem::DenseMatrix dshape;
    mfem::Array<int> dof_ids;

    for (int i = 0; i < mesh_[idx]->GetNBE(); i++) {

        if (mesh_[idx]->GetBdrAttribute(i) != face_attr) continue;

        mfem::FaceElementTransformations *FTr = mesh_[idx]->GetBdrFaceTransformations(i);
        if (FTr == nullptr) continue;

        const mfem::FiniteElement &elem = *HCurlFESpace_[idx]->GetFE(FTr->Elem1No);
        MFEM_ASSERT(elem.GetMapType() == mfem::FiniteElement::VALUE, "");
        const int int_order = 2*elem.GetOrder() + 3;
        const mfem::IntegrationRule &ir = mfem::IntRules.Get(FTr->FaceGeom, int_order);

        HCurlFESpace_[idx]->GetElementDofs(FTr->Elem1No, dof_ids);
        J_[idx]->GetSubVector(dof_ids, local_dofs);
        const int space_dim = FTr->Face->GetSpaceDim();
        normal_vec.SetSize(space_dim);
        dshape.SetSize(elem.GetDof(), space_dim);

        for (int j = 0; j < ir.GetNPoints(); j++) {

            const mfem::IntegrationPoint &ip = ir.IntPoint(j);
            mfem::IntegrationPoint eip;
            FTr->Loc1.Transform(ip, eip);
            FTr->Face->SetIntPoint(&ip);
            double face_weight = FTr->Face->Weight();
            double val = 0.0;
            FTr->Elem1->SetIntPoint(&eip);
            elem.CalcVShape(*FTr->Elem1, dshape);
            mfem::CalcOrtho(FTr->Face->Jacobian(), normal_vec);
            val += dshape.InnerProduct(normal_vec, local_dofs) / face_weight;

            // Measure the area of the boundary
            area += ip.weight * face_weight;

            // Integrate alpha * n.Grad(x) + beta * x
            flux += val * ip.weight * face_weight;
        }
      
    }

    std::cout << "Flux integration area = " << area << std::endl;
    std::cout << "Obtained flux = " << flux << std::endl;

    return flux;
}


}; // namespace hephaestus