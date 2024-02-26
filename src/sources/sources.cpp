#include "sources.hpp"

namespace hephaestus
{

void
Sources::Init(hephaestus::GridFunctions & gridfunctions,
              const hephaestus::FESpaces & fespaces,
              hephaestus::BCMap & bc_map,
              hephaestus::Coefficients & coefficients)
{
  for (const auto & [name, source] : *this)
  {
    source->Init(gridfunctions, fespaces, bc_map, coefficients);
  }
}

void
Sources::Apply(mfem::ParLinearForm * lf)
{
  for (const auto & [name, source] : *this)
  {
    source->Apply(lf);
  }
}

void
Sources::SubtractSources(mfem::ParGridFunction * gf)
{
  for (const auto & [name, source] : *this)
  {
    source->SubtractSource(gf);
  }
}

} // namespace hephaestus
