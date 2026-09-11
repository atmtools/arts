#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <format>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <variant>

namespace stdr = std::ranges;

class Wsv;

// Included by auto_wsg.h after WorkspaceGroup and WorkspaceGroupInfo are defined.
template <typename... Ts>
concept GenericContainsOnlyWorkspaceGroups = (WorkspaceGroup<std::remove_const_t<Ts>> and ...);

template <typename... Ts>
concept SortedWorkspaceGroups =
    GenericContainsOnlyWorkspaceGroups<Ts...> and
    stdr::is_sorted(std::array<std::size_t, sizeof...(Ts)>{WorkspaceGroupInfo<std::remove_const_t<Ts>>::index...});

template <typename... Ts>
concept UniformGenericConstness = (std::is_const_v<Ts> and ...) or ((not std::is_const_v<Ts>) and ...);

template <typename... Ts>
concept SortedWorkspaceGroupsUniformGenericConstness = UniformGenericConstness<Ts...> and SortedWorkspaceGroups<Ts...>;

template <SortedWorkspaceGroupsUniformGenericConstness... Ts> struct Generic : std::variant<std::shared_ptr<Ts>...> {
  using Base = std::variant<std::shared_ptr<Ts>...>;
  using Base::Base;

  Generic(const Generic&)            = default;
  Generic(Generic&&)                 = default;
  Generic& operator=(const Generic&) = default;
  Generic& operator=(Generic&&)      = default;

  // Share the contained workspace value, rejecting unsupported types at runtime.
  template <typename W> requires std::same_as<std::remove_cvref_t<W>, Wsv> Generic(W&& value) : Generic(from(value)) {}

  // All alternatives become read-only.
  using Const = Generic<const Ts...>;

 private:
  template <typename T> using Alternative = std::conditional_t<(std::same_as<T, Ts> or ...), T, const T>;

 public:
  template <typename T> requires(std::same_as<Alternative<T>, Ts> or ...) Generic(T& value)
      : Base(std::in_place_type<std::shared_ptr<Alternative<T>>>,
             std::shared_ptr<Alternative<T>>(&value, [](auto*) {})) {}

  // Temporary values may be borrowed only through const alternatives.
  template <typename T> requires(std::same_as<const T, Ts> or ...) Generic(const T& value)
      : Base(std::in_place_type<std::shared_ptr<const T>>, std::shared_ptr<const T>(&value, [](const T*) {})) {}

  template <typename T> requires(std::same_as<Alternative<T>, Ts> or ...)
  Generic(std::shared_ptr<T> value) : Base(std::in_place_type<std::shared_ptr<Alternative<T>>>, std::move(value)) {}

  // Widen an existing Generic without copying its pointee or changing ownership.
  template <typename... Us>
  requires(sizeof...(Ts) > 0) and ((std::constructible_from<Generic, std::shared_ptr<Us>>) and ...)
  Generic(const Generic<Us...>& value)
      : Generic(std::visit([](const auto& pointer) -> Generic { return pointer; }, value)) {}

  template <typename... Us>
  requires(sizeof...(Ts) > 0) and ((std::constructible_from<Generic, std::shared_ptr<Us>>) and ...)
  Generic(Generic<Us...>&& value)
      : Generic(std::visit([](auto&& pointer) -> Generic { return std::move(pointer); }, std::move(value))) {}

  template <typename W> static Generic from(const W& value) {
    std::optional<Generic> result;
    (void)((value.template holds<std::remove_const_t<Ts>>()
                ? (result.emplace(std::in_place_type<std::shared_ptr<Ts>>,
                                  value.template share<std::remove_const_t<Ts>>()),
                   true)
                : false) or
           ...);

    if (not result) {
      throw std::runtime_error(
          std::format("Unsupported generic type: {}.\n\nValid types are: {:qB,}",
                      value.type_name(),
                      std::array<std::string_view, sizeof...(Ts)>{WorkspaceGroupInfo<Ts>::name...}));
    }

    return std::move(*result);
  }
};
