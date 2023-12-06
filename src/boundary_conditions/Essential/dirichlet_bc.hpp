#pragma once
#include "essential_bc_base.hpp"

namespace hephaestus{

class DirichletBC : public EssentialBC{
public:
  DirichletBC(const std::string &name_,
                      mfem::Array<int> bdr_attributes_);
                      
  DirichletBC(const std::string &name_,
                      mfem::Array<int> bdr_attributes_,
                      mfem::ConstantCoefficient *coeff_,
                      mfem::ConstantCoefficient *coeff_im_ = nullptr);

  virtual void applyBC(mfem::GridFunction &gridfunc,
                       mfem::Mesh *mesh_) override;

  mfem::ConstantCoefficient *coeff;
  mfem::ConstantCoefficient *coeff_im;
};
} // namespace hephaestus