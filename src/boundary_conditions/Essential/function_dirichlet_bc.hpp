#pragma once
#include "essential_bc_base.hpp"

namespace hephaestus {

class FunctionDirichletBC : public EssentialBC {
public:
  FunctionDirichletBC(const std::string &name_,
                      mfem::Array<int> bdr_attributes_);
  FunctionDirichletBC(const std::string &name_,
                      mfem::Array<int> bdr_attributes_,
                      mfem::FunctionCoefficient *coeff_,
                      mfem::FunctionCoefficient *coeff_im_ = nullptr);

  virtual void applyBC(mfem::GridFunction &gridfunc,
                       mfem::Mesh *mesh_) override;

  mfem::FunctionCoefficient *coeff;
  mfem::FunctionCoefficient *coeff_im;
};

} // namespace hephaestus
