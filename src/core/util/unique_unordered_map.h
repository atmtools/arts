#pragma once

#include <format>
#include <unordered_map>

//! Class to help keep track of keys when building an unordered map
template <std::formattable<char> Key, typename T>
struct UniqueMap {
  std::unordered_map<Key, T> map;

  T& operator[](const Key& key) {
    if (map.contains(key))
      throw std::runtime_error(
          std::format("Duplicate key in UniqueMap: {}", key));

    return map[key];
  }

  const T& at(const Key& key) const { return map.at(key); }
  const T& back() const { return map.back(); }
  const T& front() const { return map.front(); }
  T& back() { return map.back(); }
  T& front() { return map.front(); }

  std::unordered_map<Key, T>::const_iterator begin() const {
    return map.begin();
  }
  std::unordered_map<Key, T>::const_iterator end() const { return map.end(); }
  std::unordered_map<Key, T>::iterator begin() { return map.begin(); }
  std::unordered_map<Key, T>::iterator end() { return map.end(); }
  std::unordered_map<Key, T>::iterator find(const Key& key) {
    return map.find(key);
  }
  std::unordered_map<Key, T>::const_iterator find(const Key& key) const {
    return map.find(key);
  }

  bool contains(const Key& key) const { return map.contains(key); }
  auto erase(std::unordered_map<Key, T>::const_iterator iter) {
    return map.erase(iter);
  }

  operator std::unordered_map<Key, T>&() { return map; }
  operator const std::unordered_map<Key, T>&() const { return map; }
};
