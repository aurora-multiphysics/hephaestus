#pragma once
#include "essential_bc_base.hpp"

namespace hephaestus {

class VectorFunctionDirichletBC : public EssentialBC {
public:
  VectorFunctionDirichletBC(const std::string &name_,
                            mfem::Array<int> bdr_attributes_);
  VectorFunctionDirichletBC(
      const std::string &name_, mfem::Array<int> bdr_attributes_,
      mfem::VectorFunctionCoefficient *vec_coeff_,
      mfem::VectorFunctionCoefficient *vec_coeff_im_ = nullptr);

  virtual void applyBC(mfem::GridFunction &gridfunc,
                       mfem::Mesh *mesh_) override;

  virtual void applyBC(mfem::ParComplexGridFunction &gridfunc,
                       mfem::Mesh *mesh_) override;

  mfem::VectorFunctionCoefficient *vec_coeff;
  mfem::VectorFunctionCoefficient *vec_coeff_im;
};

} // namespace hephaestus
