#include "variables.hpp"

namespace hephaestus {

Variable::Variable(const std::string &var_display_name_,
                   mfem::ParGridFunction *gf_)
    : var_display_name(var_display_name_), gf(gf_) {}
} // namespace hephaestus
