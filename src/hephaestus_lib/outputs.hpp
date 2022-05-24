#pragma once
#include <fstream>
#include <iostream>
#include <memory>

#include "mesh_extras.hpp"

namespace hephaestus {

class Outputs {
public:
  Outputs();
  Outputs(std::map<std::string, mfem::DataCollection*> data_collections_);

  std::map<std::string, mfem::DataCollection*> data_collections;
};

} // namespace hephaestus