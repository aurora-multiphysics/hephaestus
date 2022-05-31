#include "executioner.hpp"

namespace hephaestus {

Executioner::Executioner(const std::string &type_, const double dt_,
                         const double t_initial_, const double t_final_)
    : type(type_), dt(dt_), t_initial(t_initial_), t_final(t_final_) {}

} // namespace hephaestus
