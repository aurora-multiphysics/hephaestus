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
                    std::shared_ptr<mfem::Coefficient> coeff_,
                    std::shared_ptr<mfem::Coefficient> coeff_im_ = nullptr);

  void ApplyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_) override;

  std::shared_ptr<mfem::Coefficient> _coeff{nullptr};
  std::shared_ptr<mfem::Coefficient> _coeff_im{nullptr};
};

} // namespace hephaestus
