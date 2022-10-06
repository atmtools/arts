#pragma once

namespace Compare {
//! Returns a 'less than' lambda expression for use in, e.g., std::any_of 
auto lt(auto v) {
  return [v](const auto& x) { return x < v; };
}

//! Returns a 'less than or equal' lambda expression for use in, e.g., std::any_of
auto le(auto v) {
  return [v](const auto& x) { return x <= v; };
}

//! Returns a 'greater than' lambda expression for use in, e.g., std::any_of
auto gt(auto v) {
  return [v](const auto& x) { return x > v; };
}

//! Returns a 'greater than or equal' lambda expression for use in, e.g., std::any_of
auto ge(auto v) {
  return [v](const auto& x) { return x >= v; };
}

//! Returns a 'not equal to' lambda expression for use in, e.g., std::any_of
auto ne(auto v) {
  return [v](const auto& x) { return x not_eq v; };
}

//! Returns an 'equal to' lambda expression for use in, e.g., std::any_of
auto eq(auto v) {
  return [v](const auto& x) { return x == v; };
}
}  // namespace Compare
