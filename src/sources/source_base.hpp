#pragma once
#include "../common/pfem_extras.hpp"
#include "utils.hpp"
#include "coefficients.hpp"
#include "hephaestus_solvers.hpp"
#include "inputs.hpp"
#include "kernels.hpp"

namespace hephaestus
{

class Source : public hephaestus::Kernel<mfem::ParLinearForm>
{
public:
  virtual void SubtractSource(mfem::ParGridFunction * gf) = 0;
};

} // namespace hephaestus
