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

  // Default initializer.
  NamedFieldsMap() = default;

  // Destructor.
  ~NamedFieldsMap() { DeregisterAll(); }

  // Copy constructor. If we are creating a copy then we must NOT retain
  // ownership over any objects otherwise we will end-up with a double-free
  // situation.
  NamedFieldsMap(const NamedFieldsMap &other) {
    _field_map = other._field_map;
    _is_owner.clear(); // No ownership since copied. Prevents double-free.
  }

  // Copy assignment.
  // 1. destroy anything we have ownership of.
  // 2. Copy the field map.
  // 3. Do NOT copy the set to avoid double-freeing.
  NamedFieldsMap &operator=(const NamedFieldsMap &other) {
    DeregisterAll();
    _field_map = other._field_map;
    _is_owner.clear(); // No transferrence of ownership. Prevents double-free.

    return *this;
  }

  /// Register field @a field with name @a fname
  void Register(const std::string &fname, T *field, bool own_data) {
    // Deregister any existing field with that name and free memory if we have
    // ownership.
    Deregister(fname);

    _field_map[fname] = field;

    if (own_data) {
      _is_owner.insert(fname);
    }
  }

  /// Unregister association between field @a field and name @a fname
  void Deregister(const std::string &fname) {
    iterator iter = find(fname);
    Deregister(iter);
  }

  /// Clear all associations between names and fields
  void DeregisterAll() {
    for (iterator iter = begin(); iter != end(); iter++) {
      Deregister(iter);
    }
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
    std::vector<T *> values(keys.size());

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
  /// Deregister field and delete any owned memory.
  void Deregister(iterator iter) {
    if (iter == end()) // Not found.
      return;

    const auto &field_name = iter->first;

    if (IsDataOwner(field_name)) {
      delete iter->second;
      _is_owner.erase(field_name);
    }

    iter->second = nullptr;

    GetMap().erase(iter);
  }

  /// Returns reference to set.
  inline std::set<std::string> &GetIsDataOwnerSet() { return _is_owner; };

  /// Returns const-reference to set.
  inline const std::set<std::string> &GetIsDataOwnerSet() const {
    return _is_owner;
  };

  /// Returns true if owner of the data.
  inline bool IsDataOwner(const std::string &fname) const {
    return (GetIsDataOwnerSet().count(fname) > 0);
  }

private:
  MapType _field_map{};
  std::set<std::string> _is_owner{};
};
}; // namespace hephaestus