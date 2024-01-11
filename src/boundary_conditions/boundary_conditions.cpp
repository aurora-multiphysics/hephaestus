#include "boundary_conditions.hpp"

namespace hephaestus {

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
      if (bc != nullptr) {
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
      if (bc != nullptr) {
        bc->applyBC(gridfunc, mesh_);
      }
    }
  }
  mfem::Array<int> ess_bdr = getEssentialBdrMarkers(name_, mesh_);
  gridfunc.FESpace()->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
};

void BCMap::applyEssentialBCs(const std::string &name_,
                              mfem::Array<int> &ess_tdof_list,
                              mfem::ParComplexGridFunction &gridfunc,
                              mfem::Mesh *mesh_) {

  for (auto const &[name, bc_] : *this) {
    if (bc_->name == name_) {
      hephaestus::EssentialBC *bc =
          dynamic_cast<hephaestus::EssentialBC *>(bc_);
      if (bc != nullptr) {
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
      if (bc != nullptr) {
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
      if (bc != nullptr) {
        bc->getMarkers(*mesh_);
        bc->applyBC(clf);
      }
    }
  }
};

void BCMap::applyIntegratedBCs(const std::string &name_,
                               mfem::ParSesquilinearForm &slf,
                               mfem::Mesh *mesh_) {

  for (auto const &[name, bc_] : *this) {
    if (bc_->name == name_) {
      hephaestus::RobinBC *bc = dynamic_cast<hephaestus::RobinBC *>(bc_);
      if (bc != nullptr) {
        bc->getMarkers(*mesh_);
        bc->applyBC(slf);
      }
    }
  }
};

} // namespace hephaestus
