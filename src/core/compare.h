#pragma once

namespace Cmp {
//! Returns a 'less than' lambda expression for use in, e.g., std::any_of 
constexpr auto lt(auto v) {
  return [v](const auto& x) { return x < v; };
}

//! Returns a 'less than or equal' lambda expression for use in, e.g., std::any_of
constexpr auto le(auto v) {
  return [v](const auto& x) { return x <= v; };
}

//! Returns a 'greater than' lambda expression for use in, e.g., std::any_of
constexpr auto gt(auto v) {
  return [v](const auto& x) { return x > v; };
}

//! Returns a 'greater than or equal' lambda expression for use in, e.g., std::any_of
constexpr auto ge(auto v) {
  return [v](const auto& x) { return x >= v; };
}

//! Returns a 'not equal to' lambda expression for use in, e.g., std::any_of
constexpr auto ne(auto v) {
  return [v](const auto& x) { return x not_eq v; };
}

//! Returns an 'equal to' lambda expression for use in, e.g., std::any_of
constexpr auto eq(auto v) {
  return [v](const auto& x) { return x == v; };
}
}  // namespace Cmp
