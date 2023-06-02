#include "executioner.hpp"

namespace hephaestus {

Executioner::Executioner(const hephaestus::InputParameters &params)
    : visualization(params.GetOptionalParam<bool>("UseGLVis", false)) {}

} // namespace hephaestus
