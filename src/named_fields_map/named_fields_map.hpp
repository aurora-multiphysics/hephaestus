#pragma once
#include <map>
#include <set>
#include <string>

namespace hephaestus {

/// Lightweight adaptor over an std::map from strings to pointer to T
template <typename T> class NamedFieldsMap {
public:
  typedef std::map<std::string, T *> MapType;
  typedef typename MapType::iterator iterator;
  typedef typename MapType::const_iterator const_iterator;

  /// Default initializer.
  NamedFieldsMap() = default;

  /// Destructor.
  ~NamedFieldsMap() { DeregisterAll(); }

  /// Register association between field @a field and name @a fname
  void Register(const std::string &fname, T *field, bool own_data) {
    // 1. Deregister existing field with that name.
    Deregister(fname);

    // 2. Add to field map.
    _field_map[fname] = field;

    // 3. Keep track of ownership.
    if (own_data) {
      _owned_ptrs_map.emplace(fname, std::shared_ptr<T>(field));
    }
  }

  /// Unregister association between field @a field and name @a fname
  void Deregister(const std::string &fname) {
    iterator iter = find(fname);
    Deregister(iter);
  }

  /// Predicate to check if a field is associated with name @a fname
  inline bool Has(const std::string &fname) const {
    return find(fname) != end();
  }

  /// Get a pointer to the field associated with name @a fname
  inline T *Get(const std::string &fname) const {
    const_iterator it = find(fname);
    return it != _field_map.end() ? it->second : nullptr;
  }

  /// Returns a vector containing all values.
  std::vector<T *> Get(const std::vector<std::string> keys) {
    std::vector<T *> values;

    for (const auto &key : keys) {
      if (Has(key)) {
        values.push_back(Get(key));
      } else {
        MFEM_ABORT("Key " << key << " not found in NamedFieldsMap.");
      }
    }

    return values;
  }

  /// Returns reference to field map.
  inline MapType &GetMap() { return _field_map; }

  /// Returns const-reference to field map.
  inline const MapType &GetMap() const { return _field_map; }

  /// Returns a begin iterator to the registered fields.
  inline iterator begin() { return GetMap().begin(); }

  /// Returns a begin const iterator to the registered fields.
  inline const_iterator begin() const { return GetMap().begin(); }

  /// Returns an end iterator to the registered fields.
  inline iterator end() { return GetMap().end(); }

  /// Returns an end const iterator to the registered fields.
  inline const_iterator end() const { return GetMap().end(); }

  /// Returns an iterator to the field @a fname
  inline iterator find(const std::string &fname) {
    return GetMap().find(fname);
  }

  /// Returns a const iterator to the field @a fname
  inline const_iterator find(const std::string &fname) const {
    return GetMap().find(fname);
  }

  /// Returns the number of registered fields
  inline int NumFields() const { return GetMap().size(); }

protected:
  /// Clear all associations between names and fields
  void DeregisterAll() {
    for (iterator iter = begin(); iter != end(); iter++) {
      Deregister(iter);
    }
  }

  /// Deregister field and delete any owned memory.
  void Deregister(iterator iter) {
    if (iter == end()) // Not found.
      return;

    const auto &field_name = iter->first;

    std::cout << " Now deregistering " << field_name << std::endl;

    if (OwnsPointer(field_name)) {
      _owned_ptrs_map.erase(field_name);
    }

    _field_map.erase(field_name);
  }

  /// Returns true if we are responsible for deleting the pointer.
  bool OwnsPointer(const std::string &fname) const {
    auto iter = _owned_ptrs_map.find(fname);

    return (iter != _owned_ptrs_map.end());
  }

private:
  MapType _field_map{};
  std::map<std::string, std::shared_ptr<T>> _owned_ptrs_map{};
};
}; // namespace hephaestus