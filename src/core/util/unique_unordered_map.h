#pragma once

#include <format>
#include <unordered_map>

//! Class to help keep track of keys when building an unordered map
template <std::formattable<char> Key, typename T>
struct UniqueMap {
  std::unordered_map<Key, T> map;

  T& operator[](const Key& key) {
    if (map.contains(key)){
      throw std::runtime_error(
          std::format("Duplicate key in UniqueMap: {}", key));}

    return map[key];
  }

  const T& at(const Key& key) const { return map.at(key); }

  template <typename Self>
  decltype(auto) begin(this Self&& self) {
    return std::forward<Self>(self).map.begin();
  }

  template <typename Self>
  decltype(auto) end(this Self&& self) {
    return std::forward<Self>(self).map.end();
  }

  template <typename Self>
  decltype(auto) cbegin(this Self&& self) {
    return std::forward<Self>(self).map.cbegin();
  }

  template <typename Self>
  decltype(auto) cend(this Self&& self) {
    return std::forward<Self>(self).map.cend();
  }

  template <typename Self>
  decltype(auto) find(this Self&& self, const Key& key) {
    return std::forward<Self>(self).map.find(key);
  }

  template <typename Self>
  decltype(auto) contains(this Self&& self, const Key& key) {
    return std::forward<Self>(self).map.contains(key);
  }

  template <typename Self>
  decltype(auto) erase(this Self&& self, const Key& key) {
    return std::forward<Self>(self).map.erase(key);
  }

  template <typename Self>
  decltype(auto) erase(this Self&& self, std::unordered_map<Key, T>::const_iterator iter) {
    return std::forward<Self>(self).map.erase(iter);
  }

  template <typename Self>
  decltype(auto) erase(this Self&& self, std::unordered_map<Key, T>::iterator iter) {
    return std::forward<Self>(self).map.erase(iter);
  }

  operator std::unordered_map<Key, T>&() { return map; }
  operator const std::unordered_map<Key, T>&() const { return map; }
};
