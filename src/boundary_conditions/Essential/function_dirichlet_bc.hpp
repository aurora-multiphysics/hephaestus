#pragma once
#include "essential_bc_base.hpp"

namespace hephaestus
{

class ScalarDirichletBC : public EssentialBC
{
public:
  ScalarDirichletBC(const std::string & name_, mfem::Array<int> bdr_attributes_);
  ScalarDirichletBC(const std::string & name_,
                    mfem::Array<int> bdr_attributes_,
                    mfem::Coefficient * coeff_,
                    mfem::Coefficient * coeff_im_ = nullptr);

  void applyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_) override;

  mfem::Coefficient * coeff;
  mfem::Coefficient * coeff_im;
};

} // namespace hephaestus
