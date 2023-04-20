#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "mesh_extras.hpp"

namespace hephaestus {

class BoundaryCondition {
public:
  BoundaryCondition(const std::string &name_, mfem::Array<int> bdr_attributes_);
  mfem::Array<int> getMarkers(mfem::Mesh &mesh);

  std::string name;
  mfem::Array<int> bdr_attributes;
  mfem::Array<int> markers;

  virtual void applyBC(mfem::LinearForm &b){};
  virtual void applyBC(mfem::ComplexLinearForm &b){};
  virtual void applyBC(mfem::ParComplexLinearForm &b){};
};

class EssentialBC : public BoundaryCondition {
public:
  EssentialBC(const std::string &name_, mfem::Array<int> bdr_attributes_);

  virtual void applyBC(mfem::GridFunction &gridfunc, mfem::Mesh *mesh_){};
  virtual void applyBC(mfem::ParComplexGridFunction &gridfunc,
                       mfem::Mesh *mesh_){};
};

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

class IntegratedBC : public BoundaryCondition {
public:
  IntegratedBC(const std::string &name_, mfem::Array<int> bdr_attributes_);
  IntegratedBC(const std::string &name_, mfem::Array<int> bdr_attributes_,
               mfem::LinearFormIntegrator *lfi_re_,
               mfem::LinearFormIntegrator *lfi_im_ = nullptr);

  mfem::LinearFormIntegrator *lfi_re;
  mfem::LinearFormIntegrator *lfi_im;

  virtual void applyBC(mfem::LinearForm &b) override;
  virtual void applyBC(mfem::ComplexLinearForm &b) override;
  virtual void applyBC(mfem::ParComplexLinearForm &b) override;
};

class RobinBC : public IntegratedBC {
public:
  RobinBC(const std::string &name_, mfem::Array<int> bdr_attributes_,
          mfem::BilinearFormIntegrator *blfi_re_,
          mfem::LinearFormIntegrator *lfi_re_,
          mfem::BilinearFormIntegrator *blfi_im_ = nullptr,
          mfem::LinearFormIntegrator *lfi_im_ = nullptr);

  mfem::BilinearFormIntegrator *blfi_re;
  mfem::BilinearFormIntegrator *blfi_im;

  virtual void applyBC(mfem::ParBilinearForm &a);
  virtual void applyBC(mfem::ParSesquilinearForm &a);
};

class RWTE10PortRBC : public RobinBC {
  inline static const double epsilon0_ = 8.8541878176e-12;
  // Permeability of Free Space (units H/m)
  inline static const double mu0_ = 4.0e-7 * M_PI;

  inline static const std::complex<double> zi{std::complex<double>(0., 1.)};

public:
  RWTE10PortRBC(const std::string &name_, mfem::Array<int> bdr_attributes_,
                double frequency, double port_length_vector[3],
                double port_width_vector[3], bool input_port);

  static mfem::Vector cross_product(mfem::Vector &va, mfem::Vector &vb) {
    mfem::Vector Vec;
    Vec.SetSize(3);
    Vec[0] = va[1] * vb[2] - va[2] * vb[1];
    Vec[1] = va[2] * vb[0] - va[0] * vb[2];
    Vec[2] = va[0] * vb[1] - va[1] * vb[0];
    return Vec;
  }
  bool input_port;
  double omega_;
  mfem::Vector a1Vec;
  mfem::Vector a2Vec;
  mfem::Vector a3Vec;
  mfem::Vector a2xa3;
  mfem::Vector a3xa1;

  double V;
  double kc;
  double k0;
  std::complex<double> k_;

  mfem::Vector k_a;
  mfem::Vector k_c;

  mfem::ConstantCoefficient *robin_coef_im;
  mfem::VectorFunctionCoefficient *u_real;
  mfem::VectorFunctionCoefficient *u_imag;

  void RWTE10(const mfem::Vector &x, std::vector<std::complex<double>> &E) {

    mfem::Vector E_hat(cross_product(k_c, k_a));
    E_hat *= 1.0 / E_hat.Norml2();

    double E0(sqrt(2 * omega_ * mu0_ /
                   (a1Vec.Norml2() * a2Vec.Norml2() * k_.imag())));
    std::complex<double> E_mag =
        E0 * sin(InnerProduct(k_a, x)) * exp(-zi * InnerProduct(k_c, x));

    E[0] = E_mag * E_hat(1);
    E[1] = E_mag * E_hat(2);
    E[2] = E_mag * E_hat(0);
  }

  void RWTE10_real(const mfem::Vector &x, mfem::Vector &v) {
    std::vector<std::complex<double>> Eval(x.Size());
    RWTE10(x, Eval);
    for (int i = 0; i < x.Size(); ++i) {
      v(i) = -2 * k_.imag() * Eval[i].imag() / mu0_;
    }
  }
  void RWTE10_imag(const mfem::Vector &x, mfem::Vector &v) {
    std::vector<std::complex<double>> Eval(x.Size());
    RWTE10(x, Eval);
    for (int i = 0; i < x.Size(); ++i) {
      v(i) = 2 * k_.imag() * Eval[i].real() / mu0_;
    }
  }
};

class BCMap : public std::map<std::string, hephaestus::BoundaryCondition *> {
public:
  mfem::Array<int> getEssentialBdrMarkers(const std::string &name_,
                                          mfem::Mesh *mesh_);

  void applyEssentialBCs(const std::string &name_,
                         mfem::Array<int> &ess_tdof_list,
                         mfem::GridFunction &gridfunc, mfem::Mesh *mesh_);

  void applyEssentialBCs(const std::string &name_,
                         mfem::Array<int> &ess_tdof_list,
                         mfem::ParComplexGridFunction &gridfunc,
                         mfem::Mesh *mesh_);

  void applyIntegratedBCs(const std::string &name_, mfem::LinearForm &lf,
                          mfem::Mesh *mesh_);

  void applyIntegratedBCs(const std::string &name_,
                          mfem::ParComplexLinearForm &clf, mfem::Mesh *mesh_);

  void applyIntegratedBCs(const std::string &name_,
                          mfem::ParSesquilinearForm &clf, mfem::Mesh *mesh_);
};

} // namespace hephaestus
