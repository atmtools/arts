#pragma once

#include <concepts>
#include <type_traits>

namespace matpack {
template <typename T>
concept has_nelem = requires(T a) {
  { a.nelem() } -> std::integral;
};

template <typename T>
concept has_ncols = requires(T a) {
  { a.ncols() } -> std::integral;
};
template <typename T>
concept has_nrows = requires(T a) {
  { a.nrows() } -> std::integral;
};

template <typename T>
concept has_npages = requires(T a) {
  { a.npages() } -> std::integral;
};

template <typename T>
concept has_nbooks = requires(T a) {
  { a.nbooks() } -> std::integral;
};

template <typename T>
concept has_nshelves = requires(T a) {
  { a.nshelves() } -> std::integral;
};

template <typename T>
concept has_nvitrines = requires(T a) {
  { a.nvitrines() } -> std::integral;
};

template <typename T>
concept has_nlibraries = requires(T a) {
  { a.nlibraries() } -> std::integral;
};
}  // namespace matpack
