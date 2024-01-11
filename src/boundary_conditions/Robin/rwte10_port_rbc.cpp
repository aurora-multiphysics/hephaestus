#include "rwte10_port_rbc.hpp"

namespace hephaestus {

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

  robin_coef_im = std::make_unique<mfem::ConstantCoefficient>(k_.imag() / mu0_);
  blfi_im = std::make_unique<mfem::VectorFEMassIntegrator>(robin_coef_im.get());

  if (input_port) {
    u_real = std::make_unique<mfem::VectorFunctionCoefficient>(
        3, [this](const mfem::Vector &x, mfem::Vector &v) {
          return RWTE10_real(x, v);
        });

    u_imag = std::make_unique<mfem::VectorFunctionCoefficient>(
        3, [this](const mfem::Vector &x, mfem::Vector &v) {
          return RWTE10_imag(x, v);
        });

    lfi_re =
        std::make_unique<mfem::VectorFEBoundaryTangentLFIntegrator>(*u_real);
    lfi_im =
        std::make_unique<mfem::VectorFEBoundaryTangentLFIntegrator>(*u_imag);
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

} // namespace hephaestus
