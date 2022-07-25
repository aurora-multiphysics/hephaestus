#pragma once
#include "mfem.hpp"

namespace hephaestus {

class Variable {
public:
  Variable() {}

  Variable(const std::string &var_display_name_, mfem::ParGridFunction *gf_);

  std::string var_display_name;
  mfem::ParGridFunction *gf;
};

class VariableMap : public mfem::NamedFieldsMap<mfem::ParGridFunction> {};
} // namespace hephaestus
