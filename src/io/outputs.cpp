#include "outputs.hpp"

namespace hephaestus {

Outputs::Outputs() {}

Outputs::Outputs(
    std::map<std::string, mfem::DataCollection *> data_collections_)
    : data_collections(data_collections_) {}

} // namespace hephaestus
