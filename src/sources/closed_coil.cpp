#include "closed_coil.hpp"

namespace hephaestus {

// Pushes an element into a vector if the vector does not yet contain that same element
template<typename T> 
void pushIfUnique(std::vector<T> &vec, const T el){

    bool verify = true;

    for (auto e:vec) {
        if (e == el) verify = false;
    }

    if (verify == true) vec.push_back(el);

}

// Deletes and clears a vector of pointers
template <typename T>
void deleteAndClear(std::vector<T*> v){

    for (auto p:v) delete p;
    v.clear();
}

double highV(const mfem::Vector &x, double t) { return 1.0; }
double lowV(const mfem::Vector &x, double t) { return 0.0; }

ClosedCoilSolver::ClosedCoilSolver(const hephaestus::InputParameters &params, 
                                   const std::vector<hephaestus::Subdomain> &coil_dom,
                                   const double Jtotal,
                                   const int electrode_face,
                                   const int order) 
                                   : hcurl_fespace_name(params.GetParam<std::string>("HCurlFESpaceName")),
                                     J_gf_name(params.GetParam<std::string>("JGridFunctionName")),
                                     coil_domains_(coil_dom),
                                     Jtotal_(Jtotal),
                                     order_(order),
                                     coef1_(nullptr),
                                     coef0_(nullptr),
                                     mesh_parent_(nullptr),
                                     J_parent_(nullptr),
                                     HCurlFESpace_parent_(nullptr) {

    elec_attrs_.first = electrode_face;

}

ClosedCoilSolver::~ClosedCoilSolver() {

    delete coef1_;
    delete coef0_;
    delete high_src_;
    delete low_src_;

    deleteAndClear(mesh_);
    deleteAndClear(H1FESpace_);
    deleteAndClear(HCurlFESpace_);
    deleteAndClear(H1_Collection_);
    deleteAndClear(HCurl_Collection_);
    deleteAndClear(V_);
    deleteAndClear(sps_);
    deleteAndClear(sps_params_);
    deleteAndClear(current_solver_options_);
    deleteAndClear(gridfunctions_);
    deleteAndClear(fespaces_);
    deleteAndClear(bc_maps_);
    deleteAndClear(coefs_);
    deleteAndClear(high_DBC_);
    deleteAndClear(low_DBC_);
}

void ClosedCoilSolver::Init(hephaestus::GridFunctions &gridfunctions,
                           const hephaestus::FESpaces &fespaces, 
                           hephaestus::BCMap &bc_map,
                           hephaestus::Coefficients &coefficients) {

    coef1_ = new mfem::ConstantCoefficient(1.0);
    coef0_ = new mfem::ConstantCoefficient(0.0);

    // Retrieving the parent FE space and mesh
    HCurlFESpace_parent_ = fespaces.Get(hcurl_fespace_name);
    if (HCurlFESpace_parent_ == NULL) {
        const std::string error_message = hcurl_fespace_name +
                                        " not found in fespaces when "
                                        "creating ClosedCoilSolver\n";
        mfem::mfem_error(error_message.c_str());
    }

    J_parent_ = gridfunctions.Get(J_gf_name);
    if (J_parent_ == NULL) {
        const std::string error_message = J_gf_name +
                                        " not found in gridfunctions when "
                                        "creating ClosedCoilSolver\n";
        mfem::mfem_error(error_message.c_str());
    }


    mesh_parent_ = HCurlFESpace_parent_->GetParMesh();

    resizeChildVectors();
    makeWedge();
    initChildMeshes();
    makeFESpaces();
    makeGridFunctions();
    SetBCs();
    SPSCurrent();

    // Normalise the current through the wedges and use them as a reference
    // The MPI allreduce is to take into account the MPI partitioning
    for (int i=0; i<2; ++i) {

        double flux = calcJFlux(elec_attrs_.first, i);
        double total_flux;
        MPI_Allreduce(&flux, &total_flux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        *J_[i] /= abs(total_flux);
    }
        
}

void ClosedCoilSolver::Apply(mfem::ParLinearForm *lf) {

    for (int i=0; i<2; ++i) {
        mesh_[i]->Transfer(*J_[i], *J_parent_); 
    }

    *J_parent_ *= Jtotal_;

   lf->Add(1.0,*J_parent_);
}

void ClosedCoilSolver::SubtractSource(mfem::ParGridFunction *gf) {}

void ClosedCoilSolver::resizeChildVectors(){

    H1_Collection_.resize(2);
    HCurl_Collection_.resize(2);
    H1FESpace_.resize(2);
    HCurlFESpace_.resize(2);

    J_.resize(2);
    V_.resize(2);

    mesh_.resize(2);

    sps_params_.resize(2);
    current_solver_options_.resize(2);
    sps_.resize(2);
    gridfunctions_.resize(2);
    fespaces_.resize(2);
    bc_maps_.resize(2);
    coefs_.resize(2);

    high_DBC_.resize(2);
    low_DBC_.resize(2);

}

void ClosedCoilSolver::makeWedge(){

    std::vector<int> bdr_els;
    new_domain_attr_ = mesh_parent_->attributes.Max()+1;;
    elec_attrs_.second = mesh_parent_->bdr_attributes.Max()+1;

    // First we need to find the electrode boundary
    for (int i = 0; i<mesh_parent_->GetNBE(); ++i){
        if (mesh_parent_->GetBdrAttribute(i) == elec_attrs_.first){
            bdr_els.push_back(i);
        }
    }

    Plane3D plane;

    if (bdr_els.size() > 0){
        plane.make3DPlane(mesh_parent_, mesh_parent_->GetBdrFace(bdr_els[0]));
    }
    
    std::vector<int> elec_vtx;
    // Create a vector containing all of the vertices on the electrode
    for (auto b_fc:bdr_els) {

        mfem::Array<int> face_vtx;
        mesh_parent_->GetFaceVertices(mesh_parent_->GetBdrFace(b_fc),face_vtx);

        for (auto v:face_vtx) pushIfUnique(elec_vtx, v);
    }


    // Now we need to find all elements in the mesh that touch, on at least one vertex, the electrode face
    // if they do touch the vertex, are on one side of the electrode, and belong to the coil domain, 
    // we add them to our wedge

    std::vector<int> wedge_els;

    for (int e=0; e<mesh_parent_->GetNE(); ++e) {

        if (!isInDomain(e,coil_domains_,mesh_parent_) || 
            plane.side(elementCentre(e,mesh_parent_)) == 1) continue;

        mfem::Array<int> elem_vtx;
        mesh_parent_->GetElementVertices(e,elem_vtx);

        for (auto v1:elem_vtx) {
            for (auto v2:elec_vtx){
                if (v1 == v2) {
                        pushIfUnique(wedge_els,e);
                }
            }
        }
    }

    // Now we set the second electrode boundary attribute. Start with a list of all the 
    // faces of the wedge elements and eliminate mesh and coil boundaries, the first electrode,
    // and faces between wedge elements

    std::vector<int> wedge_faces;
    mfem::Array<int> el_faces;
    mfem::Array<int> ori;

    for (auto e:wedge_els){
        mesh_parent_->GetElementFaces(e,el_faces,ori);
        for (auto f:el_faces) pushIfUnique(wedge_faces, f);
    }

    for (auto wf:wedge_faces) {

        int e1,e2;
        mesh_parent_->GetFaceElements(wf, &e1, &e2);

        // If the face is a coil boundary
        if (!(isInDomain(e1, coil_domains_, mesh_parent_) &&
              isInDomain(e2, coil_domains_, mesh_parent_))){
            continue;
        }

        // If the face is not true interior
        if (!(mesh_parent_->FaceIsInterior(wf) || 
             (mesh_parent_->GetFaceInformation(wf).tag == mfem::Mesh::FaceInfoTag::SharedConforming ||
              mesh_parent_->GetFaceInformation(wf).tag == mfem::Mesh::FaceInfoTag::SharedSlaveNonconforming))){
            continue;
        }
        
        // If the face is shared between two elements internal to the wedge
        bool test1 = false;
        bool test2 = false;
        for (auto e:wedge_els) {
            if (e == e1) test1 = true;
            if (e == e2) test2 = true;
        }

        if (test1 && test2) continue;

        // If the face is part of the first electrode
        test1 = false;
        for (auto b_fc:bdr_els) {
            if (wf == mesh_parent_->GetBdrFace(b_fc)) {
                test1 = true;
                break;
            }
        }
        if (test1) continue;

        // At last, if the face is none of these things, it must be our second electrode
        auto* new_elem = mesh_parent_->GetFace(wf)->Duplicate(mesh_parent_);
        new_elem->SetAttribute(elec_attrs_.second);
        mesh_parent_->AddBdrElement(new_elem);
    }

    // Only after this we set the domain attributes
    for (auto e:wedge_els) mesh_parent_->SetAttribute(e, new_domain_attr_);

    mesh_parent_->FinalizeTopology();
    mesh_parent_->Finalize();
    mesh_parent_->SetAttributes();
}

void ClosedCoilSolver::initChildMeshes() {

    mfem::Array<int> doms_array;
    SubdomainToArray(coil_domains_, doms_array);
    mesh_[0] = new mfem::ParSubMesh(mfem::ParSubMesh::CreateFromDomain(*mesh_parent_,doms_array));

    hephaestus::Subdomain wedge("wedge", new_domain_attr_);
    SubdomainToArray(wedge, doms_array);
    mesh_[1] = new mfem::ParSubMesh(mfem::ParSubMesh::CreateFromDomain(*mesh_parent_,doms_array));    
    
    for (int i=0; i<2; ++i){
        inheritBdrAttributes(mesh_parent_, mesh_[i]);
    }

}

void ClosedCoilSolver::makeFESpaces(){

    for (int i=0; i<2; ++i){
        H1_Collection_[i] = new mfem::H1_FECollection(order_, mesh_[i]->Dimension());
        HCurl_Collection_[i] = new mfem::ND_FECollection(order_, mesh_[i]->Dimension());
        H1FESpace_[i] = new mfem::ParFiniteElementSpace(mesh_[i], H1_Collection_[i]);
        HCurlFESpace_[i] = new mfem::ParFiniteElementSpace(mesh_[i], HCurl_Collection_[i]);
    }
}

void ClosedCoilSolver::makeGridFunctions(){

    for (int i=0; i<2; ++i){
        V_[i] = new mfem::ParGridFunction(H1FESpace_[i]);
        J_[i] = new mfem::ParGridFunction(HCurlFESpace_[i]);
        *V_[i] = 0.0;
        *J_[i] = 0.0;
    }
}

void ClosedCoilSolver::SetBCs(){

    high_terminal_.Append(elec_attrs_.first);
    low_terminal_.Append(elec_attrs_.second);
    high_src_ = new mfem::FunctionCoefficient(highV);
    low_src_ = new mfem::FunctionCoefficient(lowV);

    for (int i=0; i<2; ++i){

        high_DBC_[i] = new hephaestus::FunctionDirichletBC(std::string("V_" + std::to_string(i)), high_terminal_, high_src_);
        low_DBC_[i] = new hephaestus::FunctionDirichletBC(std::string("V_" + std::to_string(i)), low_terminal_, low_src_);
        bc_maps_[i] = new hephaestus::BCMap;
        bc_maps_[i]->Register("high_potential", high_DBC_[i], true);
        bc_maps_[i]->Register("low_potential", low_DBC_[i], true);
    }

}

void ClosedCoilSolver::SPSCurrent(){

    for (int i=0; i<2; ++i){

        fespaces_[i] = new hephaestus::FESpaces;
        fespaces_[i]->Register(std::string("HCurl_" + std::to_string(i)), HCurlFESpace_[i], false);
        fespaces_[i]->Register(std::string("H1_" + std::to_string(i)), H1FESpace_[i], false);

        gridfunctions_[i] = new hephaestus::GridFunctions;
        gridfunctions_[i]->Register(std::string("source_" + std::to_string(i)), J_[i], false);
        gridfunctions_[i]->Register(std::string("V_" + std::to_string(i)), V_[i], false);

        current_solver_options_[i] = new hephaestus::InputParameters;
        current_solver_options_[i]->SetParam("Tolerance", float(1.0e-9));
        current_solver_options_[i]->SetParam("MaxIter", (unsigned int)1000);
        current_solver_options_[i]->SetParam("PrintLevel", 1);

        sps_params_[i] = new hephaestus::InputParameters;
        sps_params_[i]->SetParam("SourceName", std::string("source_" + std::to_string(i)));
        sps_params_[i]->SetParam("PotentialName", std::string("V_" + std::to_string(i)));
        sps_params_[i]->SetParam("HCurlFESpaceName", std::string("HCurl_" + std::to_string(i)));
        sps_params_[i]->SetParam("H1FESpaceName", std::string("H1_" + std::to_string(i)));
        sps_params_[i]->SetParam("SolverOptions", *current_solver_options_[i]);
        sps_params_[i]->SetParam("ConductivityCoefName", std::string("magnetic_permeability"));

        coefs_[i] = new hephaestus::Coefficients;
        coefs_[i]->scalars.Register("magnetic_permeability", coef1_, false);

        sps_[i] = new hephaestus::ScalarPotentialSource(*sps_params_[i]);
        sps_[i]->Init(*gridfunctions_[i], *fespaces_[i], *bc_maps_[i], *coefs_[i]);
        mfem::ParLinearForm dummy(HCurlFESpace_[i]);

        sps_[i]->Apply(&dummy);
    }

    *J_[0] *= -1.0;
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

    return flux;
}

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

bool ClosedCoilSolver::isInDomain(const int el, const std::vector<hephaestus::Subdomain> &dom, const mfem::ParMesh* mesh) {

    // This is for ghost elements
    if (el < 0) return false;

    bool verify = false;

    for (auto sd:dom) {
        if (mesh->GetAttribute(el) == sd.id) verify = true;
    }

    return verify;
}

bool ClosedCoilSolver::isInDomain(const int el, const hephaestus::Subdomain &sd, const mfem::ParMesh* mesh) {

    // This is for ghost elements
    if (el < 0) return false;

    return mesh->GetAttribute(el) == sd.id;
}

mfem::Vector ClosedCoilSolver::elementCentre(int el, mfem::ParMesh* pm){

    mfem::Array<int> elem_vtx;
    mfem::Vector com(3);
    com = 0.0;

    pm->GetElementVertices(el,elem_vtx);

    for (auto vtx:elem_vtx){
        for (int j=0; j<3; ++j) com[j] += pm->GetVertex(vtx)[j]/(double)elem_vtx.Size();
    }

    return com;
}

// 3D Plane constructor and methods

Plane3D::Plane3D():d(0) {

    u = new mfem::Vector(3);
    *u = 0.0;
}

void Plane3D::make3DPlane(const mfem::ParMesh* pm, const int face) {

    MFEM_ASSERT(pm->Dimension() == 3, "Plane3D only works in 3-dimensional meshes!");

    mfem::Array<int> face_vtx;
    std::vector<mfem::Vector> v;
    pm->GetFaceVertices(face,face_vtx);

    // First we get the coordinates of 3 vertices on the face
    for (auto vtx:face_vtx){
        mfem::Vector vtx_coords(3);
        for (int j=0; j<3; ++j) vtx_coords[j] = pm->GetVertex(vtx)[j];
        v.push_back(vtx_coords);
    }

    // Now we find the unit vector normal to the face
    v[0] -= v[1];
    v[1] -= v[2];
    v[0].cross3D(v[1],*u);
    *u /= u->Norml2();

    // Finally, we find d:
    d = *u*v[2];    
}

int Plane3D::side(const mfem::Vector v){
    double val = *u*v - d;

    if      (val > 0) return 1;
    else if (val < 0) return -1;
    else              return 0;
}


}; // namespace hephaestus