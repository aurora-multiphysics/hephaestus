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

EssentialBC::EssentialBC() {}

EssentialBC::EssentialBC(const std::string &name_,
                         mfem::Array<int> bdr_attributes_)
    : BoundaryCondition(name_, bdr_attributes_) {}

FunctionDirichletBC::FunctionDirichletBC() {}

FunctionDirichletBC::FunctionDirichletBC(const std::string &name_,
                                         mfem::Array<int> bdr_attributes_)
    : EssentialBC(name_, bdr_attributes_) {}

FunctionDirichletBC::FunctionDirichletBC(const std::string &name_,
                                         mfem::Array<int> bdr_attributes_,
                                         mfem::FunctionCoefficient *coeff_,
                                         mfem::FunctionCoefficient *coeff_im_)
    : EssentialBC(name_, bdr_attributes_), coeff(coeff_), coeff_im(coeff_im_) {}

void FunctionDirichletBC::applyBC(mfem::GridFunction &gridfunc,
                                  mfem::Mesh *mesh_) {
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = this->getMarkers(*mesh_);
  gridfunc.ProjectBdrCoefficient(*(this->coeff), ess_bdrs);
}

VectorFunctionDirichletBC::VectorFunctionDirichletBC() {}

VectorFunctionDirichletBC::VectorFunctionDirichletBC(
    const std::string &name_, mfem::Array<int> bdr_attributes_)
    : EssentialBC(name_, bdr_attributes_) {}

VectorFunctionDirichletBC::VectorFunctionDirichletBC(
    const std::string &name_, mfem::Array<int> bdr_attributes_,
    mfem::VectorFunctionCoefficient *vec_coeff_,
    mfem::VectorFunctionCoefficient *vec_coeff_im_)
    : EssentialBC(name_, bdr_attributes_), vec_coeff(vec_coeff_),
      vec_coeff_im(vec_coeff_im_) {}

void VectorFunctionDirichletBC::applyBC(mfem::GridFunction &gridfunc,
                                        mfem::Mesh *mesh_) {
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = this->getMarkers(*mesh_);
  gridfunc.ProjectBdrCoefficientTangent(*(this->vec_coeff), ess_bdrs);
}

IntegratedBC::IntegratedBC() {}

IntegratedBC::IntegratedBC(const std::string &name_,
                           mfem::Array<int> bdr_attributes_)
    : BoundaryCondition(name_, bdr_attributes_) {}

IntegratedBC::IntegratedBC(const std::string &name_,
                           mfem::Array<int> bdr_attributes_,
                           mfem::LinearFormIntegrator *lfi_re_,
                           mfem::LinearFormIntegrator *lfi_im_)
    : BoundaryCondition(name_, bdr_attributes_), lfi_re(lfi_re_),
      lfi_im(lfi_im_) {}

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
  global_ess_markers = 0;
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = 0;
  hephaestus::EssentialBC *bc;
  for (auto const &[name, bc_] : *this) {
    if (bc_->name == name_) {
      bc = dynamic_cast<hephaestus::EssentialBC *>(bc_);
      if (bc != NULL) {
        ess_bdrs = bc->getMarkers(*mesh_);
        for (auto it = 0; it != mesh_->bdr_attributes.Max(); ++it) {
          global_ess_markers[it] =
              std::max(global_ess_markers[it], ess_bdrs[it]);
        }
      }
    }
  }
  return global_ess_markers;
}

void BCMap::applyEssentialBCs(const std::string &name_,
                              mfem::Array<int> &ess_tdof_list,
                              mfem::GridFunction &gridfunc, mfem::Mesh *mesh_) {

  for (auto const &[name, bc_] : *this) {
    if (bc_->name == name_) {
      hephaestus::EssentialBC *bc =
          dynamic_cast<hephaestus::EssentialBC *>(bc_);
      if (bc != NULL) {
        bc->applyBC(gridfunc, mesh_);
      }
    }
  }
  mfem::Array<int> ess_bdr = getEssentialBdrMarkers(name_, mesh_);
  gridfunc.FESpace()->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
};

void BCMap::applyIntegratedBCs(const std::string &name_, mfem::LinearForm &lf,
                               mfem::Mesh *mesh_) {

  for (auto const &[name, bc_] : *this) {
    if (bc_->name == name_) {
      hephaestus::IntegratedBC *bc =
          dynamic_cast<hephaestus::IntegratedBC *>(bc_);
      if (bc != NULL) {
        bc->getMarkers(*mesh_);
        bc->applyBC(lf);
      }
    }
  }
};

void BCMap::applyIntegratedBCs(const std::string &name_,
                               mfem::ParComplexLinearForm &clf,
                               mfem::Mesh *mesh_) {

  for (auto const &[name, bc_] : *this) {
    if (bc_->name == name_) {
      hephaestus::IntegratedBC *bc =
          dynamic_cast<hephaestus::IntegratedBC *>(bc_);
      if (bc != NULL) {
        bc->getMarkers(*mesh_);
        bc->applyBC(clf);
      }
    }
  }
};
} // namespace hephaestus
