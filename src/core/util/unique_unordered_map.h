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

  operator std::unordered_map<Key, T>() && { return std::move(map); }
};
