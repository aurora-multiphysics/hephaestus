#include "rwte10_port_rbc.hpp"

namespace hephaestus
{

RWTE10PortRBC::RWTE10PortRBC(const std::string & name_,
                             mfem::Array<int> bdr_attributes_,
                             double frequency_,
                             double port_length_vector_[3],
                             double port_width_vector_[3],
                             bool input_port_)
  : RobinBC(name_, bdr_attributes_, nullptr, nullptr, nullptr, nullptr),
    input_port(input_port_),
    omega_(2 * M_PI * frequency_),
    a1Vec(port_length_vector_, 3),
    a2Vec(port_width_vector_, 3),
    a3Vec(CrossProduct(a1Vec, a2Vec)),
    a2xa3(CrossProduct(a2Vec, a3Vec)),
    a3xa1(CrossProduct(a3Vec, a1Vec)),
    V(mfem::InnerProduct(a1Vec, a2xa3)),
    kc(M_PI / a1Vec.Norml2()),
    k0(omega_ * sqrt(epsilon0_ * mu0_)),
    k_(std::complex<double>(0., sqrt(k0 * k0 - kc * kc))),
    k_a(a2xa3),
    k_c(a3Vec)
{
  k_a *= M_PI / V;
  k_c *= k_.imag() / a3Vec.Norml2();

  robin_coef_im = std::make_unique<mfem::ConstantCoefficient>(k_.imag() / mu0_);
  blfi_im = std::make_unique<mfem::VectorFEMassIntegrator>(robin_coef_im.get());

  if (input_port)
  {
    u_real = std::make_unique<mfem::VectorFunctionCoefficient>(
        3, [this](const mfem::Vector & x, mfem::Vector & v) { return RWTE10Real(x, v); });

    u_imag = std::make_unique<mfem::VectorFunctionCoefficient>(
        3, [this](const mfem::Vector & x, mfem::Vector & v) { return RWTE10Imag(x, v); });

    lfi_re = std::make_unique<mfem::VectorFEBoundaryTangentLFIntegrator>(*u_real);
    lfi_im = std::make_unique<mfem::VectorFEBoundaryTangentLFIntegrator>(*u_imag);
  }
}

void
RWTE10PortRBC::RWTE10(const mfem::Vector & x, std::vector<std::complex<double>> & E)
{

  mfem::Vector e_hat(CrossProduct(k_c, k_a));
  e_hat *= 1.0 / e_hat.Norml2();

  double e0(sqrt(2 * omega_ * mu0_ / (a1Vec.Norml2() * a2Vec.Norml2() * k_.imag())));
  std::complex<double> e_mag = e0 * sin(InnerProduct(k_a, x)) * exp(-zi * InnerProduct(k_c, x));

  E[0] = e_mag * e_hat(1);
  E[1] = e_mag * e_hat(2);
  E[2] = e_mag * e_hat(0);
}

void
RWTE10PortRBC::RWTE10Real(const mfem::Vector & x, mfem::Vector & v)
{
  std::vector<std::complex<double>> eval(x.Size());
  RWTE10(x, eval);
  for (int i = 0; i < x.Size(); ++i)
  {
    v(i) = -2 * k_.imag() * eval[i].imag() / mu0_;
  }
}
void
RWTE10PortRBC::RWTE10Imag(const mfem::Vector & x, mfem::Vector & v)
{
  std::vector<std::complex<double>> eval(x.Size());
  RWTE10(x, eval);
  for (int i = 0; i < x.Size(); ++i)
  {
    v(i) = 2 * k_.imag() * eval[i].real() / mu0_;
  }
}

} // namespace hephaestus
