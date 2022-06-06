#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "mesh_extras.hpp"

namespace hephaestus {

class BoundaryCondition {
public:
  BoundaryCondition();
  BoundaryCondition(const std::string &name_, mfem::Array<int> bdr_attributes_);
  mfem::Array<int> getMarkers(mfem::Mesh &mesh);

  std::string name;
  mfem::Array<int> bdr_attributes;
  mfem::Array<int> markers;

  virtual void applyBC(mfem::LinearForm &b){};
  virtual void applyBC(mfem::ComplexLinearForm &b){};
  virtual void applyBC(mfem::ParComplexLinearForm &b){};
};

class FunctionDirichletBC : public BoundaryCondition {
public:
  FunctionDirichletBC();
  FunctionDirichletBC(const std::string &name_,
                      mfem::Array<int> bdr_attributes_);
  FunctionDirichletBC(const std::string &name_,
                      mfem::Array<int> bdr_attributes_,
                      mfem::FunctionCoefficient *coeff_,
                      mfem::FunctionCoefficient *coeff_im_ = nullptr);

  mfem::FunctionCoefficient *coeff;
  mfem::FunctionCoefficient *coeff_im;
};

class VectorFunctionDirichletBC : public BoundaryCondition {
public:
  VectorFunctionDirichletBC();
  VectorFunctionDirichletBC(const std::string &name_,
                            mfem::Array<int> bdr_attributes_);
  VectorFunctionDirichletBC(
      const std::string &name_, mfem::Array<int> bdr_attributes_,
      mfem::VectorFunctionCoefficient *vec_coeff_,
      mfem::VectorFunctionCoefficient *vec_coeff_im_ = nullptr);

  mfem::VectorFunctionCoefficient *vec_coeff;
  mfem::VectorFunctionCoefficient *vec_coeff_im;
};

class IntegratedBC : public BoundaryCondition {
public:
  IntegratedBC();
  IntegratedBC(const std::string &name_, mfem::Array<int> bdr_attributes_);

  mfem::LinearFormIntegrator *lfi_re;
  mfem::LinearFormIntegrator *lfi_im;

  virtual void applyBC(mfem::LinearForm &b) override;
  virtual void applyBC(mfem::ComplexLinearForm &b) override;
  virtual void applyBC(mfem::ParComplexLinearForm &b) override;
};

class BCMap : public std::map<std::string, hephaestus::BoundaryCondition *> {
public:
  //  global_ess_markers;
  mfem::Array<int> getEssentialBdrMarkers(const std::string &name_,
                                          mfem::Mesh *mesh_) {

    mfem::Array<int> global_ess_markers(mesh_->bdr_attributes.Max());
    mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
    hephaestus::FunctionDirichletBC *bc;
    for (auto const &[name, bc_] : *this) {
      if (bc_->name == name_) {
        bc = dynamic_cast<hephaestus::FunctionDirichletBC *>(bc_);
        ess_bdrs = bc->getMarkers(*mesh_);
        for (auto it = 0; it != mesh_->bdr_attributes.Max(); ++it) {
          global_ess_markers[it] =
              std::max(global_ess_markers[it], ess_bdrs[it]);
        }
      }
      // bc->coeff->SetTime(this->GetTime());
      // Phi_gf.ProjectBdrCoefficient(*(bc->coeff), ess_bdrs);
    }
    return global_ess_markers;
  };

}; // namespace hephaestus

} // namespace hephaestus
