#pragma once

#include <cstddef>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>

#include "debug.h"

namespace FieldMap {
template <typename T>
concept Hashable = requires(T a) {
  {
    std::hash<std::remove_cvref_t<T>>{}(a)
  } -> std::convertible_to<std::size_t>;
};

template <typename T>
concept Equatable = requires(T a) {
  { a == a } -> std::convertible_to<bool>;
};

template <typename T>
concept field_key = Hashable<T> and Equatable<T>;

template <typename T, typename... Ts>
concept one_of = (std::same_as<std::remove_cvref_t<T>, Ts> or ...);

/** A multi-key map of fields of data
 *
 * This will create something that should behave similar to
 * a std::unordered_map<std::variant<Keys...>, T>.  If something
 * can be done by std::unordered_map but not by this type, please
 * feel free to add it.
 * 
 * @tparam T The type of data
 * @tparam Keys Any number of keys so that std::unordered_map<Key, T> is possible
 */
template <typename T, field_key... Keys>
  requires(sizeof...(Keys) > 0)
struct Map {
  static_assert(not(std::same_as<Keys, bool> or ...),
                "Cannot have pure boolean keys");
  static_assert((std::same_as<Keys, std::remove_cvref_t<Keys>> and ...),
                "Only for base-type keys");

  static constexpr std::size_t N = sizeof...(Keys);
  using KeyVal = std::variant<Keys...>;
  std::tuple<std::unordered_map<Keys, T>...> map_data;

 private:
  template <one_of<Keys...> Key>
  static constexpr std::size_t pos() {
    constexpr std::array<bool, N> tst{
        std::same_as<std::remove_cvref_t<Key>, Keys>...};
    std::size_t i = 0;
    while (i < N and not tst[i]) i++;
    return i;
  }

 public:
  template <one_of<Keys...> Key>
  std::unordered_map<std::remove_cvref_t<Key>, T> &map()
    requires(pos<Key>() < N)
  {
    return std::get<pos<Key>()>(map_data);
  }

  template <one_of<Keys...> Key>
  const std::unordered_map<std::remove_cvref_t<Key>, T> &map() const
    requires(pos<Key>() < N)
  {
    return std::get<pos<Key>()>(map_data);
  }

  template <one_of<Keys...> Key>
  [[nodiscard]] constexpr const T &operator[](const Key &k) const try {
    return map<Key>().at(k);
  } catch (std::out_of_range &) {
    throw std::out_of_range(
        var_string("Key not found in map: ", std::quoted(var_string(k))));
  } catch (...) {
    throw;
  }

  template <one_of<Keys...> Key>
  [[nodiscard]] constexpr T &operator[](const Key &k) {
    return map<Key>()[k];
  }

  [[nodiscard]] constexpr const T &operator[](const KeyVal &k) const try {
    return std::visit(
        [this](auto &key) -> const T & {
          return this->map<decltype(key)>().at(key);
        },
        k);
  } catch (std::out_of_range &) {
    throw std::out_of_range(
        var_string("Key not found in map: ", std::quoted(var_string(k))));
  } catch (...) {
    throw;
  }

  [[nodiscard]] constexpr T &operator[](const KeyVal &k) try {
    return std::visit(
        [this](auto &key) -> T & {
          return const_cast<Map *>(this)->map<decltype(key)>()[key];
        },
        k);
  }
  ARTS_METHOD_ERROR_CATCH

  template <one_of<Keys...> Key>
  constexpr void erase(Key &&k) {
    map<Key>().erase(std::forward<Key>(k));
  }

  template <one_of<Keys...> Key>
  constexpr bool contains(Key &&key) const {
    return map<Key>().contains(std::forward<Key>(key));
  }

  constexpr bool contains(const KeyVal &key) const {
    return std::visit(
        [this](auto &k) -> bool {
          return this->map<decltype(k)>().contains(k);
        },
        key);
  }

  template <field_key... MultipleKeys>
  constexpr bool has(MultipleKeys &&...key) const {
    return (contains<MultipleKeys>(std::forward<MultipleKeys>(key)) and ...);
  }

 private:
  template <one_of<Keys...> Key, field_key... rest_of_keys>
  constexpr void append_keys(auto &keys) const {
    for (auto &val : map<Key>()) keys.emplace_back(val.first);
    if constexpr (sizeof...(rest_of_keys) > 0)
      append_keys<rest_of_keys...>(keys);
  }

 public:
  template <one_of<Keys..., bool> Key = bool>
  constexpr auto keys() const try {
    if constexpr (std::same_as<Key, bool>) {
      std::vector<KeyVal> out;
      out.reserve(size());
      append_keys<Keys...>(out);
      return out;
    } else {
      std::vector<Key> out;
      out.reserve(size<Key>());
      append_keys<Key>(out);
      return out;
    }
  }
  ARTS_METHOD_ERROR_CATCH

  template <one_of<Keys..., bool> Key = bool>
  [[nodiscard]] constexpr std::size_t size() const {
    if constexpr (std::same_as<Key, bool>) {
      return (size<Keys>() + ...);
    } else {
      return map<Key>().size();
    }
  }

  template <one_of<Keys...> Key>
  [[nodiscard]] auto begin() {
    return map<Key>().begin();
  }

  template <one_of<Keys...> Key>
  [[nodiscard]] auto end() {
    return map<Key>().end();
  }

  template <one_of<Keys...> Key>
  [[nodiscard]] auto begin() const {
    return map<Key>().begin();
  }

  template <one_of<Keys...> Key>
  [[nodiscard]] auto end() const {
    return map<Key>().end();
  }

  template <one_of<Keys...> Key>
  [[nodiscard]] auto cbegin() const {
    return map<Key>().begin();
  }

  template <one_of<Keys...> Key>
  [[nodiscard]] auto cend() const {
    return map<Key>().end();
  }
};
}  // namespace FieldMap
