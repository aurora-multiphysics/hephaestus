#pragma once
#include <map>
#include <set>
#include <string>

namespace hephaestus
{

/// Lightweight adaptor over an std::map from strings to pointer to T
template <typename T>
class NamedFieldsMap
{
public:
  typedef std::map<std::string, T *> MapType;
  typedef typename MapType::iterator iterator;
  typedef typename MapType::const_iterator const_iterator;

  /// Default initializer.
  NamedFieldsMap() = default;

  /// Destructor.
  ~NamedFieldsMap() { DeregisterAll(); }

  /// Register association between field @a field and name @a field_name
  void Register(const std::string & field_name, T * field, bool own_data)
  {
    // 1. Deregister existing field with that name.
    Deregister(field_name);

    // 2. Add to field map.
    _field_map[field_name] = field;

    // 3. Keep track of ownership.
    if (own_data)
    {
      _stored_ptrs_map.emplace(field_name, std::shared_ptr<T>(field));
    }
  }

  /// Unregister association between field @a field and name @a field_name.
  void Deregister(const std::string & field_name)
  {
    _field_map.erase(field_name);
    _stored_ptrs_map.erase(field_name);
  }

  /// Predicate to check if a field is associated with name @a field_name.
  inline bool Has(const std::string & field_name) const { return find(field_name) != end(); }

  /// Get a pointer to the field associated with name @a field_name.
  inline T * Get(const std::string & field_name) const
  {
    const_iterator it = find(field_name);
    return it != _field_map.end() ? it->second : nullptr;
  }

  /// Returns a vector containing all values for supplied keys.
  std::vector<T *> Get(const std::vector<std::string> keys)
  {
    std::vector<T *> values;

    for (const auto & key : keys)
    {
      if (Has(key))
        values.push_back(Get(key));
      else
        MFEM_ABORT("Key " << key << " not found in NamedFieldsMap.");
    }

    values.shrink_to_fit();
    return values;
  }

  /// Returns reference to field map.
  inline MapType & GetMap() { return _field_map; }

  /// Returns const-reference to field map.
  inline const MapType & GetMap() const { return _field_map; }

  /// Returns a begin iterator to the registered fields.
  inline iterator begin() { return _field_map.begin(); }

  /// Returns a begin const iterator to the registered fields.
  inline const_iterator begin() const { return _field_map.begin(); }

  /// Returns an end iterator to the registered fields.
  inline iterator end() { return _field_map.end(); }

  /// Returns an end const iterator to the registered fields.
  inline const_iterator end() const { return _field_map.end(); }

  /// Returns an iterator to the field @a field_name.
  inline iterator find(const std::string & field_name) { return _field_map.find(field_name); }

  /// Returns a const iterator to the field @a field_name.
  inline const_iterator find(const std::string & field_name) const
  {
    return _field_map.find(field_name);
  }

  /// Returns the number of registered fields.
  inline int NumFields() const { return _field_map.size(); }

protected:
  /// Clear all associations between names and fields.
  void DeregisterAll()
  {
    _field_map.clear();
    _stored_ptrs_map.clear();
  }

private:
  MapType _field_map{};
  std::map<std::string, std::shared_ptr<T>> _stored_ptrs_map{};
};
}; // namespace hephaestus