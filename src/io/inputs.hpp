#pragma once
#include <any>
#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_map>
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
  InputParameters(std::unordered_map<std::string, std::any> params) : _params(std::move(params)) {}

  /// Removes all parameters.
  void Clear() { _params.clear(); }

  /// Sets an existing parameter.
  template <typename T>
  void Set(const std::string & name, T & value)
  {
    SetParam(name, std::any(value));
  }

  /// Sets an existing parameter.
  void SetParam(const std::string & name, std::any value) { _params[name] = value; }

  /// Returns a copy of the stored parameter.
  template <typename T>
  [[nodiscard]] T Get(const std::string & name) const
  {
    CheckForMissingParam(name);

    return std::any_cast<T>(_params.at(name));
  }

  /// Returns the stored parameter or a default if not found.
  template <typename T>
  [[nodiscard]] T GetOptionalParam(const std::string & name, T default_value) const
  {
    return Has(name) ? Get<T>(name) : default_value;
  }

  /// Returns true if there exists a parameter matching the name.
  [[nodiscard]] bool Has(const std::string & name) const
  {
    auto it = _params.find(name);
    return (it != _params.end());
  }

protected:
  /// Checks for a missing parameter.
  void CheckForMissingParam(const std::string & name) const
  {
    if (!Has(name))
    {
      MFEM_ABORT("No parameter with name '" << name << "'.");
    }
  }

private:
  std::unordered_map<std::string, std::any> _params;
};

} // namespace hephaestus
