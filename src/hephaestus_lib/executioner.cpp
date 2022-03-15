#include "executioner.hpp"

namespace hephaestus {

Executioner::Executioner(const std::string& executioner_type,
                         const double time_step, const double end_time)
    : type(executioner_type), dt(time_step), t_final(end_time) {}

}  // namespace hephaestus