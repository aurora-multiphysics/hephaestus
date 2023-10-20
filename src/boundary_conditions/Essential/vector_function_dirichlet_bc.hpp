#pragma once
#include "essential_bc_base.hpp"

namespace hephaestus {

class VectorFunctionDirichletBC : public EssentialBC {
public:
  VectorFunctionDirichletBC(const std::string &name_,
                            mfem::Array<int> bdr_attributes_);
  VectorFunctionDirichletBC(
      const std::string &name_, mfem::Array<int> bdr_attributes_,
      mfem::VectorCoefficient *vec_coeff_,
      mfem::VectorCoefficient *vec_coeff_im_ = nullptr);

  virtual void applyBC(mfem::GridFunction &gridfunc,
                       mfem::Mesh *mesh_) override;

  virtual void applyBC(mfem::ParComplexGridFunction &gridfunc,
                       mfem::Mesh *mesh_) override;

  mfem::VectorCoefficient *vec_coeff;
  mfem::VectorCoefficient *vec_coeff_im;
};

} // namespace hephaestus
