#pragma once

#include <concepts>
#include <format>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <variant>

template <typename... Ts>
concept UniformGenericConstness = (std::is_const_v<Ts> and ...) or ((not std::is_const_v<Ts>) and ...);

template <typename... Ts> requires UniformGenericConstness<Ts...> struct Generic;

namespace workspace_variant_detail {
template <typename Result, typename... Others> struct Extend;
template <typename... Ts> struct Extend<Generic<Ts...>> {
  using type = Generic<Ts...>;
};
template <typename... Ts, typename... Rest> struct Extend<Generic<Ts...>, Generic<>, Rest...>
    : Extend<Generic<Ts...>, Rest...> {};
template <typename... Ts, typename U, typename... Us, typename... Rest>
struct Extend<Generic<Ts...>, Generic<U, Us...>, Rest...>
    : Extend<std::conditional_t<(std::same_as<U, Ts> or ...), Generic<Ts...>, Generic<Ts..., U>>,
             Generic<Us...>,
             Rest...> {};
}  // namespace workspace_variant_detail

// Alternatives retain workspace ownership. Constness belongs to each pointee.
// Direct C++ references are borrowed and must not escape the call.
template <typename... Ts> requires UniformGenericConstness<Ts...>
struct Generic : std::variant<std::shared_ptr<Ts>...> {
  using Base = std::variant<std::shared_ptr<Ts>...>;
  using Base::Base;
  // Ordered union of alternative types; exact duplicates are removed.
  template <typename... Others> using Extend = typename workspace_variant_detail::Extend<Generic, Others...>::type;

  // All alternatives become read-only.
  using Const = typename workspace_variant_detail::Extend<Generic<>, Generic<const Ts>...>::type;

  [[nodiscard]] Const as_const() const& { return Const(*this); }
  [[nodiscard]] Const as_const() && { return Const(std::move(*this)); }

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

  // Lvalue variants borrow their active value. Do not destroy or change the
  // source alternative while the Generic is in use.
  template <typename... Us> requires(sizeof...(Ts) > 0) and ((std::constructible_from<Generic, Us&>) and ...)
  Generic(std::variant<Us...>& value) : Generic(std::visit([](auto& item) -> Generic { return item; }, value)) {}

  template <typename... Us> requires(sizeof...(Ts) > 0) and ((std::constructible_from<Generic, const Us&>) and ...)
  Generic(const std::variant<Us...>& value)
      : Generic(std::visit([](const auto& item) -> Generic { return item; }, value)) {}

  // Own a moved value so a temporary variant cannot leave a dangling pointer.
  template <typename... Us>
  requires(sizeof...(Ts) > 0) and ((std::constructible_from<Generic, std::shared_ptr<Us>>) and ...)
  Generic(std::variant<Us...>&& value)
      : Generic(std::visit(
            [](auto&& item) -> Generic {
              using T = std::remove_cvref_t<decltype(item)>;
              return std::make_shared<T>(std::move(item));
            },
            std::move(value))) {}

  template <typename... Us>
  requires(sizeof...(Ts) > 0) and ((std::constructible_from<Generic, std::shared_ptr<const Us>>) and ...)
  Generic(const std::variant<Us...>&& value)
      : Generic(std::visit(
            [](const auto& item) -> Generic {
              using T = std::remove_cvref_t<decltype(item)>;
              return std::make_shared<const T>(item);
            },
            value)) {}

  template <typename... Us> Generic(std::variant<Us...>&&)       = delete;
  template <typename... Us> Generic(const std::variant<Us...>&&) = delete;

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
    if (not result) throw std::runtime_error(std::format("Unsupported generic type: {}", value.type_name()));
    return std::move(*result);
  }
};

// Deduction mirrors value-variant constness; rvalues transfer their active value.
template <typename... Ts> Generic(std::variant<Ts...>&) -> Generic<Ts...>;
template <typename... Ts> Generic(const std::variant<Ts...>&) -> Generic<const Ts...>;
template <typename... Ts> Generic(std::variant<Ts...>&&) -> Generic<Ts...>;
template <typename... Ts> Generic(const std::variant<Ts...>&&) -> Generic<const Ts...>;
