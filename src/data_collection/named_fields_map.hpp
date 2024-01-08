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
  ~NamedFieldsMap() { DeleteData(); }

  /// Register field @a field with name @a fname
  void Register(const std::string &fname, T *field, bool own_data) {
    T *&ref = _field_map[fname];

    if (own_data) {
      delete ref; // if newly allocated -> ref is null -> OK
      GetIsDataOwnerSet().insert(fname);
    }

    ref = field;
  }

  /// Unregister association between field @a field and name @a fname
  void Deregister(const std::string &fname) {
    iterator iter = find(fname);
    Deregister(iter);
  }

  /// Clear all associations between names and fields
  void DeleteData() {
    for (iterator iter = begin(); iter != end(); iter++) {
      DeleteData(iter);
    }

    GetMap().clear();
  }

  /// Predicate to check if a field is associated with name @a fname
  inline bool Has(const std::string &fname) const {
    return find(fname) != end();
  }

  /// Get a pointer to the field associated with name @a fname
  inline T *Get(const std::string &fname) const {
    const_iterator it = find(fname);
    return it != _field_map.end() ? it->second : nullptr;
  };

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

    DeleteData(iter);
    GetMap().erase(iter);
  }

  /// Deletes memory for field.
  void DeleteData(iterator iter) {
    const auto &field_name = iter->first;

    if (IsDataOwner(field_name)) {
      delete iter->second;
      GetIsDataOwnerSet().erase(field_name);
    }

    iter->second = nullptr;
  }

  /// Returns reference to set.
  inline std::set<std::string> &GetIsDataOwnerSet() {
    return _is_data_owner_set;
  };

  /// Returns const-reference to set.
  inline const std::set<std::string> &GetIsDataOwnerSet() const {
    return _is_data_owner_set;
  };

  /// Returns true if owner of the data.
  inline bool IsDataOwner(const std::string &fname) const {
    return (GetIsDataOwnerSet().count(fname) > 0);
  }

private:
  MapType _field_map{};
  std::set<std::string> _is_data_owner_set{};
};
}; // namespace hephaestus