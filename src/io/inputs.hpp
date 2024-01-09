#pragma once
#include <any>
#include <fstream>
#include <iostream>
#include <memory>

#include "boundary_conditions.hpp"
#include "coefficients.hpp"
#include "outputs.hpp"

namespace hephaestus
{

class InputParameters
{

protected:
  std::map<std::string, std::any> params;

public:
  InputParameters(){};
  InputParameters(std::map<std::string, std::any> _params) : params(_params){};
  void SetParam(std::string param_name, std::any value) { params[param_name] = value; };
  template <typename T>
  T GetParam(std::string param_name) const
  {
    T param;
    try
    {
      param = std::any_cast<T>(params.at(param_name));
    }
    catch (const std::exception & e)
    {
      MFEM_ABORT("Exception raised when trying to cast required parameter " << param_name);
    }
    return param;
  };
  template <typename T>
  T GetOptionalParam(std::string param_name, T value) const
  {
    T param;
    try
    {
      param = std::any_cast<T>(params.at(param_name));
    }
    catch (...)
    {
      param = value;
    }
    return param;
  };
};

} // namespace hephaestus
