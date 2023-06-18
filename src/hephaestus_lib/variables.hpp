#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.hpp"
#include "mfem.hpp"

namespace hephaestus {

class FESpaces : public mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> {};

class GridFunctions : public mfem::NamedFieldsMap<mfem::ParGridFunction> {};

} // namespace hephaestus
