#pragma once
#include <any>
#include <fstream>
#include <iostream>
#include <memory>
#include <utility>

#include "boundary_conditions.hpp"
#include "coefficients.hpp"
#include "outputs.hpp"

namespace hephaestus
{

class InputParameters
{
public:
  InputParameters() = default;
  InputParameters(std::map<std::string, std::any> params) : _params(std::move(params)) {}

  /// Adds a new parameter matching the name.
  template <typename T>
  void AddParam(const std::string & name, T & value)
  {
    _params[name] = std::any(value);
  }

  /// Returns a reference to the parameter matching the name.
  template <typename T>
  [[nodiscard]] T & Set(const std::string & name)
  {
    return std::any_cast<T &>(_params.at(name));
  }

  /// Returns true if there exists a parameter matching the name.
  [[nodiscard]] bool Has(const std::string & name) const
  {
    auto it = _params.find(name);
    return (it != _params.end());
  }

  void SetParam(std::string param_name, std::any value) { _params[param_name] = value; }

  template <typename T>
  [[nodiscard]] T GetParam(std::string param_name) const
  {
    T param;

    try
    {
      param = std::any_cast<T>(_params.at(param_name));
    }
    catch (const std::exception & e)
    {
      MFEM_ABORT("Exception raised when trying to cast required parameter '" << param_name
                                                                             << "': " << e.what());
    }

    return param;
  }

  template <typename T>
  [[nodiscard]] T GetOptionalParam(std::string param_name, T value) const
  {
    T param;

    try
    {
      param = std::any_cast<T>(_params.at(param_name));
    }
    catch (...)
    {
      param = value;
    }

    return param;
  }

protected:
  std::map<std::string, std::any> _params;
};

} // namespace hephaestus
