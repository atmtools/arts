#pragma once

namespace Cmp {
//! Returns a 'less than' lambda expression for use in, e.g., std::any_of
constexpr auto lt(auto v) {
  return [v](const auto& x) { return x < v; };
}

template <auto V>
constexpr auto lt() {
  return [](const auto& x) static { return x < V; };
}

//! Returns a 'less than or equal' lambda expression for use in, e.g., std::any_of
constexpr auto le(auto v) {
  return [v](const auto& x) { return x <= v; };
}

template <auto V>
constexpr auto le() {
  return [](const auto& x) static { return x <= V; };
}

//! Returns a 'greater than' lambda expression for use in, e.g., std::any_of
constexpr auto gt(auto v) {
  return [v](const auto& x) { return x > v; };
}

template <auto V>
constexpr auto gt() {
  return [](const auto& x) static { return x > V; };
}

//! Returns a 'greater than or equal' lambda expression for use in, e.g., std::any_of
constexpr auto ge(auto v) {
  return [v](const auto& x) { return x >= v; };
}

template <auto V>
constexpr auto ge() {
  return [](const auto& x) static { return x >= V; };
}

//! Returns a 'not equal to' lambda expression for use in, e.g., std::any_of
constexpr auto ne(auto v) {
  return [v](const auto& x) { return x != v; };
}

template <auto V>
constexpr auto ne() {
  return [](const auto& x) static { return x != V; };
}

//! Returns an 'equal to' lambda expression for use in, e.g., std::any_of
constexpr auto eq(auto v) {
  return [v](const auto& x) { return x == v; };
}

template <auto V>
constexpr auto eq() {
  return [](const auto& x) static { return x == V; };
}

//! Returns a `this->contains` lambda expression for use in, e.g., std::any_of
constexpr auto contains(const auto& v) {
  return [&v](const auto& x) { return x.contains(v); };
}

template <auto V>
constexpr auto contains() {
  return [](const auto& x) static { return x.contains(V); };
}
}  // namespace Cmp
