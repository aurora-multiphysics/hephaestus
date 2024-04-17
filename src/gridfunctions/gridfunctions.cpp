#include "gridfunctions.hpp"

namespace hephaestus
{

void
FESpaces::Update()
{
  for ([[maybe_unused]] const auto & [name, fespace] : *this)
  {
    fespace->Update();
  }
}

void
GridFunctions::Update()
{
  for ([[maybe_unused]] const auto & [name, gridfunction] : *this)
  {
    gridfunction->Update();
  }
}

}