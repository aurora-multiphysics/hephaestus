#include "named_fields_map.hpp"

using namespace hephaestus;

template <typename T> NamedFieldsMap<T>::~NamedFieldsMap<T>() { DeleteData(); }

template <typename T>
void NamedFieldsMap<T>::Register(const std::string &fname, T *field,
                                 bool own_data) {
  T *&ref = _field_map[fname];

  if (own_data) {
    delete ref; // if newly allocated -> ref is null -> OK
    GetIsDataOwnerSet().insert(fname);
  }

  ref = field;
}

template <typename T>
void NamedFieldsMap<T>::Deregister(const std::string &fname) {
  iterator iter = find(fname);
  Deregister(iter);
}

template <typename T> void NamedFieldsMap<T>::Deregister(iterator iter) {
  if (iter == end()) // Not found.
    return;

  DeleteData(iter);
  GetMap().erase(iter);
}

template <typename T> void NamedFieldsMap<T>::DeleteData(iterator iter) {
  const auto &field_name = iter->first;

  if (IsDataOwner(field_name)) {
    delete iter->second;
    GetIsDataOwnerSet().erase(field_name);
  }

  iter->second = nullptr;
}

template <typename T> void NamedFieldsMap<T>::DeleteData() {
  for (iterator iter = begin(); iter != end(); iter++) {
    DeleteData(iter);
  }

  GetMap().clear();
}
