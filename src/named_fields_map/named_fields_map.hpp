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

  /// Construct new field with name field_name and register.
  template <class FieldType, class... FieldArgs>
  void Register(const std::string & field_name, FieldArgs &&... args)
  {
    Register(field_name, std::make_shared<FieldType>(std::forward<FieldArgs>(args)...));
  }

  /// Register association between field and field_name.
  void Register(const std::string & field_name, std::shared_ptr<T> field)
  {
    CheckForNonNullField(field_name, field.get());
    CheckForDoubleRegistration(field_name, field.get());

    Deregister(field_name);

    _field_map[field_name] = std::move(field);
  }

  /// Unregister association between a field and the field name.
  void Deregister(const std::string & field_name) { _field_map.erase(field_name); }

  /// Predicate to check if a field is associated with name field_name.
  [[nodiscard]] inline bool Has(const std::string & field_name) const
  {
    return find(field_name) != end();
  }

  /// Returns a shared pointer to the field associated with name field_name.
  [[nodiscard]] inline std::shared_ptr<T> GetShared(const std::string & field_name) const
  {
    CheckForRegistration(field_name);

    return find(field_name)->second;
  }

  /// Returns a reference to a field.
  [[nodiscard]] inline T & GetRef(const std::string & field_name) const
  {
    return *GetShared(field_name);
  }

  /// Returns a non-owning pointer to the field associated with name field_name.
  [[nodiscard]] inline T * Get(const std::string & field_name) const
  {
    return GetShared(field_name).get();
  }

  /// Returns a non-owning pointer to the field where TDerived is a derived class of class T.
  template <typename TDerived>
  [[nodiscard]] inline TDerived * Get(const std::string & field_name) const
  {
    auto ptr = Get(field_name);

    return dynamic_cast<TDerived *>(ptr);
  }

  /// Returns a vector containing all values for supplied keys.
  [[nodiscard]] std::vector<T *> Get(const std::vector<std::string> keys) const
  {
    std::vector<T *> values;

    for (const auto & key : keys)
    {
      values.push_back(Get(key));
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
    if (Has(field_name) && Get(field_name) == field)
    {
      MFEM_ABORT("The field '" << field_name << "' is already registered.");
    }
  }

  /// Check that a field with that name exists in the map.
  void CheckForRegistration(const std::string & field_name) const
  {
    if (!Has(field_name))
    {
      MFEM_ABORT("The field '" << field_name << "' has not been registered.");
    }
  }

  /// Check that the field to be registered is non-null.
  void CheckForNonNullField(const std::string & field_name, T * field) const
  {
    if (!field)
    {
      MFEM_ABORT("Cannot register NULL field with name '" << field_name << "'.");
    }
  }

  /// Clear all associations between names and fields.
  void DeregisterAll() { _field_map.clear(); }

private:
  MapType _field_map{};
};
} // namespace hephaestus