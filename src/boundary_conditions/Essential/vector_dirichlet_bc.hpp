#pragma once
#include "essential_bc_base.hpp"

namespace hephaestus
{

class VectorDirichletBC : public EssentialBC
{

protected:
  enum APPLY_TYPE
  {
    STANDARD,
    TANGENTIAL,
    NORMAL
  };

public:
  VectorDirichletBC(const std::string & name_, mfem::Array<int> bdr_attributes_);
  VectorDirichletBC(const std::string & name_,
                    mfem::Array<int> bdr_attributes_,
                    mfem::VectorCoefficient * vec_coeff_,
                    mfem::VectorCoefficient * vec_coeff_im_ = nullptr,
                    APPLY_TYPE boundary_apply_type = TANGENTIAL);

  virtual void applyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_) override;

  virtual void applyBC(mfem::ParComplexGridFunction & gridfunc, mfem::Mesh * mesh_) override;

  mfem::VectorCoefficient * vec_coeff;
  mfem::VectorCoefficient * vec_coeff_im;
  APPLY_TYPE boundary_apply_type;
};

} // namespace hephaestus
