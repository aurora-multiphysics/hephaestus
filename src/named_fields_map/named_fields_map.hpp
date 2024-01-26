#pragma once
#include <map>
#include <set>
#include <string>
#include <memory>
#include <vector>
#include <mfem.hpp>

namespace hephaestus
{

/// Lightweight adaptor over an std::map from strings to pointer to T
template <typename T>
class NamedFieldsMap
{
public:
  using MapType = std::map<std::string, std::shared_ptr<T>>;
  using iterator = typename MapType::iterator;
  using const_iterator = typename MapType::const_iterator;

  /// Default initializer.
  NamedFieldsMap() = default;

  /// Destructor.
  ~NamedFieldsMap() { DeregisterAll(); }

  /// Construct new field with name @field_name and register.
  template <class FieldType, class... FieldArgs>
  void Register(const std::string & field_name, FieldArgs &&... args)
  {
    Register(field_name, std::make_shared<FieldType>(std::forward<FieldArgs>(args)...));
  }

  /// Register association between field @a field and name @a field_name.
  void Register(const std::string & field_name, std::shared_ptr<T> field)
  {
    CheckForDoubleRegistration(field_name, field.get());

    Deregister(field_name);

    _field_map[field_name] = std::move(field);
  }

  /// Unregister association between field @a field and name @a field_name.
  void Deregister(const std::string & field_name) { _field_map.erase(field_name); }

  /// Predicate to check if a field is associated with name @a field_name.
  [[nodiscard]] inline bool Has(const std::string & field_name) const
  {
    return find(field_name) != end();
  }

  /// Get a shared pointer to the field associated with name @a field_name.
  [[nodiscard]] inline std::shared_ptr<T> Get(const std::string & field_name) const
  {
    auto it = find(field_name);
    return it != _field_map.end() ? it->second : nullptr;
  }

  /// Get a shared pointer to the field and cast to subclass TDerived.
  template <typename TDerived>
  [[nodiscard]] inline std::shared_ptr<TDerived> Get(const std::string & field_name) const
  {
    auto owned_ptr = Get(field_name);

    return owned_ptr ? std::dynamic_pointer_cast<TDerived>(owned_ptr) : nullptr;
  }

  /// Returns a vector containing all values for supplied keys.
  std::vector<std::shared_ptr<T>> Get(const std::vector<std::string> keys) const
  {
    std::vector<std::shared_ptr<T>> values;

    for (const auto & key : keys)
    {
      if (Has(key))
        values.push_back(Get(key));
      else
      {
        std::string key_not_found_msg("Key " + key + " not found in NamedFieldsMap.");
        MFEM_ABORT(key_not_found_msg);
      }
    }

    values.shrink_to_fit();
    return values;
  }

  /// Returns reference to field map.
  inline MapType & GetMap() { return _field_map; }

  /// Returns const-reference to field map.
  [[nodiscard]] inline const MapType & GetMap() const { return _field_map; }

  /// Returns a begin iterator to the registered fields.
  // NOLINTNEXTLINE(readability-identifier-naming)
  inline iterator begin() { return _field_map.begin(); }

  /// Returns a begin const iterator to the registered fields.
  // NOLINTNEXTLINE(readability-identifier-naming)
  [[nodiscard]] inline const_iterator begin() const { return _field_map.begin(); }

  /// Returns an end iterator to the registered fields.
  // NOLINTNEXTLINE(readability-identifier-naming)
  inline iterator end() { return _field_map.end(); }

  /// Returns an end const iterator to the registered fields.
  // NOLINTNEXTLINE(readability-identifier-naming)
  [[nodiscard]] inline const_iterator end() const { return _field_map.end(); }

  /// Returns an iterator to the field @a field_name.
  // NOLINTNEXTLINE(readability-identifier-naming)
  inline iterator find(const std::string & field_name) { return _field_map.find(field_name); }

  /// Returns a const iterator to the field @a field_name.
  // NOLINTNEXTLINE(readability-identifier-naming)
  [[nodiscard]] inline const_iterator find(const std::string & field_name) const
  {
    return _field_map.find(field_name);
  }

  /// Returns the number of registered fields.
  [[nodiscard]] inline int NumFields() const { return _field_map.size(); }

protected:
  /// Check for double-registration of a field. A double-registered field may
  /// result in undefined behavior.
  void CheckForDoubleRegistration(const std::string & field_name, T * field) const
  {
    if (Has(field_name) && Get(field_name).get() == field)
    {
      MFEM_ABORT("The field '" << field_name << "' is already registered.");
    }
  }

  /// Clear all associations between names and fields.
  void DeregisterAll() { _field_map.clear(); }

private:
  MapType _field_map{};
};
} // namespace hephaestus