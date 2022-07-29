#pragma once
#include <any>
#include <fstream>
#include <iostream>
#include <memory>

#include "boundary_conditions.hpp"
#include "executioner.hpp"
#include "materials.hpp"
#include "outputs.hpp"

namespace hephaestus {

class InputParameters {

protected:
  std::map<std::string, std::any> params;

public:
  InputParameters(){};
  void SetParam(std::string param_name, std::any value) {
    params[param_name] = value;
  };
  template <typename T> T GetParam(std::string param_name) const {
    T param;
    try {
      param = std::any_cast<T>(params.at(param_name));
    } catch (const std::bad_any_cast &e) {
      std::cout << e.what() << '\n';
    }
    return param;
  };
};

class Inputs {
public:
  Inputs(){};
  Inputs(const mfem::Mesh &mesh_, const std::string &formulation_,
         const int order_, const BCMap &bc_map_,
         const DomainProperties &domain_properties_,
         const Executioner &executioner_, Outputs outputs_);

  mfem::Mesh mesh;
  std::string formulation;
  int order;
  BCMap bc_map;
  DomainProperties domain_properties;
  Executioner executioner;
  Outputs outputs;
};

} // namespace hephaestus
