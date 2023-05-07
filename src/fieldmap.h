#pragma once

#include "matpack_concepts.h"

#include <concepts>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>

namespace FieldMap {
template<typename T>
concept Hashable = requires(T a) {
    { std::hash<std::remove_cvref_t<T>>{}(a) } -> std::convertible_to<std::size_t>;
};

template <typename T>
concept Equatable = requires(T a) {
  { a == a } -> std::convertible_to<bool>;
};

template <typename T> concept field_key = Hashable<T> and Equatable<T>;

/** A multi-key map of fields of data
 *
 * This will create something that should behave similar to
 * a std::unordered_map<std::variant<Keys...>, T>.  If something
 * can be done by std::unordered_map but not by this type, please
 * feel free to add it.
 *
 * In addition to that, it should interact well with ARTS concepts,
 * such as Index being the unit of a nelem.
 * 
 * @tparam T The type of data
 * @tparam Keys Any number of keys so that std::unordered_map<Key, T> is possible
 */
template <typename T, field_key... Keys>
  requires(sizeof...(Keys) > 0)
struct Map {
static_assert(not (std::same_as<Keys, bool> or ...), "Cannot have pure boolean keys");
static_assert((std::same_as<Keys, std::remove_cvref_t<Keys>> and ...), "Only for base-type keys");

  static constexpr Index N = sizeof...(Keys);
  using KeyVal = std::variant<Keys...>;
  std::tuple<std::unordered_map<Keys, T>...> map_data;

private:
  template <field_key Key> static constexpr Index pos() {
    constexpr std::array<bool, N> tst{
        std::same_as<std::remove_cvref_t<Key>, Keys>...};
    Index i=0;
    while(not tst[i] and i < N) i++;
    return i;
  }

public:
  template <field_key Key>
  std::unordered_map<std::remove_cvref_t<Key>, T> &map()
    requires(pos<Key>() < N)
  {
    return std::get<pos<Key>()>(map_data);
  }

  template <field_key Key>
  const std::unordered_map<std::remove_cvref_t<Key>, T> &map() const
    requires(pos<Key>() < N)
  {
    return std::get<pos<Key>()>(map_data);
  }

  template <field_key Key>
  [[nodiscard]] constexpr const T &operator[](const Key &k) const {
    return map<Key>().at(k);
  }

  template <field_key Key>
  [[nodiscard]] constexpr T &operator[](const Key &k) {
    return map<Key>()[k];
  }

  [[nodiscard]] constexpr const T &operator[](const KeyVal &k) const {
    return std::visit([this](auto& key) -> const T& {
      return this->map<decltype(key)>().at(key);
      }, k);
  }

  [[nodiscard]] constexpr T &operator[](const KeyVal &k) {
    return std::visit([this](auto& key) -> T& {
      return const_cast<Map *>(this)->map<decltype(key)>()[key];
      }, k);
  }

  template <field_key Key> constexpr void erase(Key &&k) {
    map<Key>().erase(std::forward<Key>(k));
  }

  template <field_key Key>
  constexpr bool contains(Key &&key) const {
    return map<Key>().contains(std::forward<Key>(key));
  }

  template <field_key ... MultipleKeys>
  constexpr bool has(MultipleKeys && ...key) const {
    return (contains<MultipleKeys>(std::forward<MultipleKeys>(key)) and ...);
  }

private:
  template <field_key Key, field_key... rest_of_keys>
  constexpr void append_keys(auto &keys) const {
    for (auto &val : map<Key>())
      keys.emplace_back(val.first);
    if constexpr (sizeof...(rest_of_keys) > 0)
      append_keys<rest_of_keys...>(keys);
  }

public:
  template <field_key Key = bool> constexpr auto keys() const {
    if constexpr (std::same_as<Key, bool>) {
      std::vector<KeyVal> out;
      append_keys<Keys...>(out);
      return out;
    } else {
      std::vector<Key> out;
      append_keys<Key>(out);
      return out;
    }
  }

  template <field_key Key = bool> [[nodiscard]] constexpr Index nelem() const {
    if constexpr (std::same_as<Key, bool>) {
      return (nelem<Keys>() + ...);
    } else {
      return static_cast<Index>(map<Key>().size());
    }
  }

  template <field_key Key>
  [[nodiscard]] auto begin() {return map<Key>().begin();}

  template <field_key Key>
  [[nodiscard]] auto end() {return map<Key>().end();}

  template <field_key Key>
  [[nodiscard]] auto begin() const {return map<Key>().begin();}

  template <field_key Key>
  [[nodiscard]] auto end() const {return map<Key>().end();}

  template <field_key Key>
  [[nodiscard]] auto cbegin() const {return map<Key>().begin();}

  template <field_key Key>
  [[nodiscard]] auto cend() const {return map<Key>().end();}
};
} // namespace FieldMap
