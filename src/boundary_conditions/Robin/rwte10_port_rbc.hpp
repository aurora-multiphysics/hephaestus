#pragma once
#include "robin_bc_base.hpp"

namespace hephaestus {

class RWTE10PortRBC : public RobinBC {
  inline static const double epsilon0_{8.8541878176e-12};
  inline static const double mu0_{4.0e-7 * M_PI};

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
  void RWTE10(const mfem::Vector &x, std::vector<std::complex<double>> &E);
  void RWTE10_real(const mfem::Vector &x, mfem::Vector &v);
  void RWTE10_imag(const mfem::Vector &x, mfem::Vector &v);

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

  std::unique_ptr<mfem::ConstantCoefficient> robin_coef_im;
  std::unique_ptr<mfem::VectorFunctionCoefficient> u_real;
  std::unique_ptr<mfem::VectorFunctionCoefficient> u_imag;
};

} // namespace hephaestus
