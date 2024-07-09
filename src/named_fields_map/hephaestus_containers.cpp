#include "hephaestus_containers.hpp"
#include "hephaestus.hpp"

namespace hephaestus
{

void
FESpaces::Update()
{
  for ([[maybe_unused]] const auto & [name, fespace] : *this)
  {
    logger.debug("Update called for fespace '{}'.", name);
    fespace->Update();
  }
}

void
GridFunctions::Update()
{
  for ([[maybe_unused]] const auto & [name, gridfunction] : *this)
  {
    logger.debug("Update called for gridfunction '{}'.", name);
    gridfunction->Update();
  }
}

}