#pragma once
#include "mfem.hpp"

namespace hephaestus {

class FECollections
    : public mfem::NamedFieldsMap<mfem::FiniteElementCollection> {};

class FESpaces : public mfem::NamedFieldsMap<mfem::ParFiniteElementSpace> {};

class GridFunctions : public mfem::NamedFieldsMap<mfem::ParGridFunction> {};

} // namespace hephaestus
