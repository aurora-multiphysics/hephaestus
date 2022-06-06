#include "boundary_conditions.hpp"

namespace hephaestus {

BoundaryCondition::BoundaryCondition() {}

BoundaryCondition::BoundaryCondition(const std::string &name_,
                                     mfem::Array<int> bdr_attributes_)
    : name(name_), bdr_attributes(bdr_attributes_) {}

mfem::Array<int> BoundaryCondition::getMarkers(mfem::Mesh &mesh) {
  mfem::common::AttrToMarker(mesh.bdr_attributes.Max(), bdr_attributes,
                             markers);
  return markers;
}

FunctionDirichletBC::FunctionDirichletBC() {}

FunctionDirichletBC::FunctionDirichletBC(const std::string &name_,
                                         mfem::Array<int> bdr_attributes_)
    : BoundaryCondition(name_, bdr_attributes_) {}

FunctionDirichletBC::FunctionDirichletBC(const std::string &name_,
                                         mfem::Array<int> bdr_attributes_,
                                         mfem::FunctionCoefficient *coeff_,
                                         mfem::FunctionCoefficient *coeff_im_)
    : BoundaryCondition(name_, bdr_attributes_), coeff(coeff_),
      coeff_im(coeff_im_) {}

VectorFunctionDirichletBC::VectorFunctionDirichletBC() {}

VectorFunctionDirichletBC::VectorFunctionDirichletBC(
    const std::string &name_, mfem::Array<int> bdr_attributes_)
    : BoundaryCondition(name_, bdr_attributes_) {}

VectorFunctionDirichletBC::VectorFunctionDirichletBC(
    const std::string &name_, mfem::Array<int> bdr_attributes_,
    mfem::VectorFunctionCoefficient *vec_coeff_,
    mfem::VectorFunctionCoefficient *vec_coeff_im_)
    : BoundaryCondition(name_, bdr_attributes_), vec_coeff(vec_coeff_),
      vec_coeff_im(vec_coeff_im_) {}

IntegratedBC::IntegratedBC() {}

IntegratedBC::IntegratedBC(const std::string &name_,
                           mfem::Array<int> bdr_attributes_)
    : BoundaryCondition(name_, bdr_attributes_) {}

void IntegratedBC::applyBC(mfem::LinearForm &b) {
  b.AddBoundaryIntegrator(lfi_re, markers);
}

void IntegratedBC::applyBC(mfem::ComplexLinearForm &b) {
  b.AddBoundaryIntegrator(lfi_re, lfi_im, markers);
}

void IntegratedBC::applyBC(mfem::ParComplexLinearForm &b) {
  b.AddBoundaryIntegrator(lfi_re, lfi_im, markers);
}

mfem::Array<int> BCMap::getEssentialBdrMarkers(const std::string &name_,
                                               mfem::Mesh *mesh_) {
  mfem::Array<int> global_ess_markers(mesh_->bdr_attributes.Max());
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  hephaestus::FunctionDirichletBC *bc;
  for (auto const &[name, bc_] : *this) {
    if (bc_->name == name_) {
      bc = dynamic_cast<hephaestus::FunctionDirichletBC *>(bc_);
      ess_bdrs = bc->getMarkers(*mesh_);
      for (auto it = 0; it != mesh_->bdr_attributes.Max(); ++it) {
        global_ess_markers[it] = std::max(global_ess_markers[it], ess_bdrs[it]);
      }
    }
  }
  return global_ess_markers;
}

mfem::Array<int> BCMap::applyEssentialBCs(const std::string &name_,
                                          mfem::GridFunction &gridfunc,
                                          mfem::Mesh *mesh_, double time) {

  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());

  for (auto const &[name, bc_] : *this) {
    if (bc_->name == name_) {
      hephaestus::FunctionDirichletBC *bc =
          dynamic_cast<hephaestus::FunctionDirichletBC *>(bc_);
      ess_bdrs = bc->getMarkers(*mesh_);
      bc->coeff->SetTime(time);
      gridfunc.ProjectBdrCoefficient(*(bc->coeff), ess_bdrs);
    }
  }
  return getEssentialBdrMarkers(name_, mesh_);
};

} // namespace hephaestus
