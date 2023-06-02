#include "boundary_conditions.hpp"

namespace hephaestus {

BoundaryCondition::BoundaryCondition(const std::string &name_,
                                     mfem::Array<int> bdr_attributes_)
    : name(name_), bdr_attributes(bdr_attributes_) {}

mfem::Array<int> BoundaryCondition::getMarkers(mfem::Mesh &mesh) {
  mfem::common::AttrToMarker(mesh.bdr_attributes.Max(), bdr_attributes,
                             markers);
  return markers;
}

EssentialBC::EssentialBC(const std::string &name_,
                         mfem::Array<int> bdr_attributes_)
    : BoundaryCondition(name_, bdr_attributes_) {}

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
  if (this->vec_coeff == NULL) {
    MFEM_ABORT(
        "Boundary condition does not store valid coefficients to specify the "
        "components of the vector at the Dirichlet boundary.");
  }
  gridfunc.ProjectBdrCoefficientTangent(*(this->vec_coeff), ess_bdrs);
}

void VectorFunctionDirichletBC::applyBC(mfem::ParComplexGridFunction &gridfunc,
                                        mfem::Mesh *mesh_) {
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = this->getMarkers(*mesh_);
  if (this->vec_coeff == NULL || this->vec_coeff_im == NULL) {
    MFEM_ABORT(
        "Boundary condition does not store valid coefficients to specify both "
        "the real and imaginary components of the vector at the Dirichlet "
        "boundary.");
  }
  gridfunc.ProjectBdrCoefficientTangent(*(this->vec_coeff),
                                        *(this->vec_coeff_im), ess_bdrs);
}

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

RobinBC::RobinBC(const std::string &name_, mfem::Array<int> bdr_attributes_,
                 mfem::BilinearFormIntegrator *blfi_re_,
                 mfem::LinearFormIntegrator *lfi_re_,
                 mfem::BilinearFormIntegrator *blfi_im_,
                 mfem::LinearFormIntegrator *lfi_im_)
    : IntegratedBC(name_, bdr_attributes_, lfi_re_, lfi_im_), blfi_re(blfi_re_),
      blfi_im(blfi_im_) {}

void RobinBC::applyBC(mfem::ParBilinearForm &a) {
  a.AddBoundaryIntegrator(blfi_re, markers);
}

void RobinBC::applyBC(mfem::ParSesquilinearForm &a) {
  a.AddBoundaryIntegrator(blfi_re, blfi_im, markers);
}

RWTE10PortRBC::RWTE10PortRBC(const std::string &name_,
                             mfem::Array<int> bdr_attributes_,
                             double frequency_, double port_length_vector_[3],
                             double port_width_vector_[3], bool input_port_)
    : RobinBC(name_, bdr_attributes_, NULL, NULL, NULL, NULL),
      input_port(input_port_), omega_(2 * M_PI * frequency_),
      a1Vec(port_length_vector_, 3), a2Vec(port_width_vector_, 3),
      a3Vec(cross_product(a1Vec, a2Vec)), a2xa3(cross_product(a2Vec, a3Vec)),
      a3xa1(cross_product(a3Vec, a1Vec)), V(mfem::InnerProduct(a1Vec, a2xa3)),
      kc(M_PI / a1Vec.Norml2()), k0(omega_ * sqrt(epsilon0_ * mu0_)),
      k_(std::complex<double>(0., sqrt(k0 * k0 - kc * kc))), k_a(a2xa3),
      k_c(a3Vec) {
  k_a *= M_PI / V;
  k_c *= k_.imag() / a3Vec.Norml2();

  robin_coef_im = new mfem::ConstantCoefficient(k_.imag() / mu0_);
  blfi_im = new mfem::VectorFEMassIntegrator(robin_coef_im);

  if (input_port) {
    u_real = new mfem::VectorFunctionCoefficient(
        3, [this](const mfem::Vector &x, mfem::Vector &v) {
          return RWTE10_real(x, v);
        });

    u_imag = new mfem::VectorFunctionCoefficient(
        3, [this](const mfem::Vector &x, mfem::Vector &v) {
          return RWTE10_imag(x, v);
        });

    lfi_re = new mfem::VectorFEBoundaryTangentLFIntegrator(*u_real);
    lfi_im = new mfem::VectorFEBoundaryTangentLFIntegrator(*u_imag);
  }
}

void RWTE10PortRBC::RWTE10(const mfem::Vector &x,
                           std::vector<std::complex<double>> &E) {

  mfem::Vector E_hat(cross_product(k_c, k_a));
  E_hat *= 1.0 / E_hat.Norml2();

  double E0(
      sqrt(2 * omega_ * mu0_ / (a1Vec.Norml2() * a2Vec.Norml2() * k_.imag())));
  std::complex<double> E_mag =
      E0 * sin(InnerProduct(k_a, x)) * exp(-zi * InnerProduct(k_c, x));

  E[0] = E_mag * E_hat(1);
  E[1] = E_mag * E_hat(2);
  E[2] = E_mag * E_hat(0);
}

void RWTE10PortRBC::RWTE10_real(const mfem::Vector &x, mfem::Vector &v) {
  std::vector<std::complex<double>> Eval(x.Size());
  RWTE10(x, Eval);
  for (int i = 0; i < x.Size(); ++i) {
    v(i) = -2 * k_.imag() * Eval[i].imag() / mu0_;
  }
}
void RWTE10PortRBC::RWTE10_imag(const mfem::Vector &x, mfem::Vector &v) {
  std::vector<std::complex<double>> Eval(x.Size());
  RWTE10(x, Eval);
  for (int i = 0; i < x.Size(); ++i) {
    v(i) = 2 * k_.imag() * Eval[i].real() / mu0_;
  }
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

void BCMap::applyEssentialBCs(const std::string &name_,
                              mfem::Array<int> &ess_tdof_list,
                              mfem::ParComplexGridFunction &gridfunc,
                              mfem::Mesh *mesh_) {

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

void BCMap::applyIntegratedBCs(const std::string &name_,
                               mfem::ParSesquilinearForm &slf,
                               mfem::Mesh *mesh_) {

  for (auto const &[name, bc_] : *this) {
    if (bc_->name == name_) {
      hephaestus::RobinBC *bc = dynamic_cast<hephaestus::RobinBC *>(bc_);
      if (bc != NULL) {
        bc->getMarkers(*mesh_);
        bc->applyBC(slf);
      }
    }
  }
};

} // namespace hephaestus
