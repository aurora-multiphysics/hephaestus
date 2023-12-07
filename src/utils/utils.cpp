#include "utils.hpp"

namespace hephaestus {

double calcFlux(mfem::GridFunction *v_field, int face_attr) {

  double flux = 0.0;
  double area = 0.0;

  mfem::FiniteElementSpace *FES = v_field->FESpace();
  mfem::Mesh *mesh = FES->GetMesh();

  mfem::Vector local_dofs, normal_vec;
  mfem::DenseMatrix dshape;
  mfem::Array<int> dof_ids;

  for (int i = 0; i < mesh->GetNBE(); i++) {

    if (mesh->GetBdrAttribute(i) != face_attr)
      continue;

    mfem::FaceElementTransformations *FTr =
        mesh->GetFaceElementTransformations(mesh->GetBdrElementFaceIndex(i));
    if (FTr == nullptr)
      continue;

    const mfem::FiniteElement &elem = *FES->GetFE(FTr->Elem1No);
    const int int_order = 2 * elem.GetOrder() + 3;
    const mfem::IntegrationRule &ir =
        mfem::IntRules.Get(FTr->FaceGeom, int_order);

    FES->GetElementDofs(FTr->Elem1No, dof_ids);
    v_field->GetSubVector(dof_ids, local_dofs);
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

  double total_flux;
  MPI_Allreduce(&flux, &total_flux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return total_flux;
}

void SubdomainToArray(const std::vector<hephaestus::Subdomain> &sd,
                      mfem::Array<int> &arr) {
  arr.DeleteAll();
  for (auto s : sd)
    arr.Append(s.id);
}

void SubdomainToArray(const hephaestus::Subdomain &sd, mfem::Array<int> &arr) {
  arr.DeleteAll();
  arr.Append(sd.id);
}

void inheritBdrAttributes(const mfem::ParMesh *parent_mesh,
                          mfem::ParSubMesh *child_mesh) {

  int face, ori, att;
  auto map = child_mesh->GetParentToSubMeshFaceIDMap();

  for (int bdr = 0; bdr < parent_mesh->GetNBE(); ++bdr) {

    parent_mesh->GetBdrElementFace(bdr, &face, &ori);
    if (map[face] != -1) {
      att = parent_mesh->GetBdrAttribute(bdr);
      auto *new_elem = child_mesh->GetFace(map[face])->Duplicate(child_mesh);
      new_elem->SetAttribute(att);
      child_mesh->AddBdrElement(new_elem);
    }
  }

  child_mesh->FinalizeTopology();
  child_mesh->Finalize();
  child_mesh->SetAttributes();
}

void attrToMarker(const mfem::Array<int> attr_list,
                  mfem::Array<int> &marker_list, int max_attr) {

  marker_list.SetSize(max_attr);
  marker_list = 0;

  for (auto a : attr_list)
    marker_list[a - 1] = 1;
}

void cleanDivergence(mfem::ParGridFunction &Vec_GF,
                     hephaestus::InputParameters solve_pars) {

  hephaestus::InputParameters pars;
  hephaestus::GridFunctions gfs;
  hephaestus::FESpaces fes;
  hephaestus::BCMap bcs;

  gfs.Register("Vector_GF", &Vec_GF, false);
  pars.SetParam("VectorGridFunctionName", std::string("Vector_GF"));
  pars.SetParam("SolverOptions", solve_pars);
  hephaestus::HelmholtzProjector projector(pars);
  projector.Project(gfs, fes, bcs);
}

void cleanDivergence(hephaestus::GridFunctions &gfs, hephaestus::BCMap &bcs,
                     const std::string vec_gf_name,
                     const std::string scalar_gf_name,
                     hephaestus::InputParameters solve_pars) {

  hephaestus::InputParameters pars;
  hephaestus::FESpaces fes;

  pars.SetParam("VectorGridFunctionName", vec_gf_name);
  pars.SetParam("ScalarGridFunctionName", scalar_gf_name);
  pars.SetParam("SolverOptions", solve_pars);
  hephaestus::HelmholtzProjector projector(pars);
  projector.Project(gfs, fes, bcs);
}

} // namespace hephaestus