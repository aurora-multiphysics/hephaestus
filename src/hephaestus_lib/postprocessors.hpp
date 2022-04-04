#include "mfem.hpp"

namespace hephaestus {

class JouleHeatingCoefficient : public mfem::Coefficient {
private:
  mfem::ParGridFunction &E_gf;
  mfem::PWCoefficient sigma;

public:
  JouleHeatingCoefficient(const mfem::PWCoefficient &sigma_,
                          mfem::ParGridFunction &E_gf_)
      : E_gf(E_gf_), sigma(sigma_) {}
  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip);
  virtual ~JouleHeatingCoefficient() {}
};

} // namespace hephaestus
