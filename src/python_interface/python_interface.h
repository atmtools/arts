#ifndef python_interface_h
#define python_interface_h

#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <variant>
#include <optional>

#include <auto_md.h>
#include <debug.h>
#include <enums.h>

//! Contains a bunch of helper functions to manipulate python objects inside C++
namespace Python {
namespace py = pybind11;

struct Numeric_ {
  Numeric val;
  constexpr operator Numeric&() { return val; }
  constexpr operator const Numeric&() const { return val; }
  constexpr Numeric_& operator=(Numeric x) {val = x; return *this;}
};

struct Index_ {
  Index val;
  constexpr operator Index&() { return val; }
  constexpr operator const Index&() const { return val; };
  constexpr Index_& operator=(Index x) {val = x; return *this;}
};

template <class T, class VariantT>
T& select_gout(const VariantT& val) {
  return std::visit([](auto&& out) -> T& {return *out;}, val);
}

template <class T, class VariantT>
const T& select_gin(const VariantT& val) {
  return std::visit([](auto&& out) -> const T& {return *out;}, val);
}

template <class T, class VariantT>
const T& select_gin(const T& default_, const std::optional<VariantT>& val) {
  if (val) return std::visit([](auto&& out) -> const T& {return *out;}, val.value());
  return default_;
}

template <class T, class WorkspaceVariable, class VariantT>
T& select_out(WorkspaceVariable wsv, std::optional<VariantT>& val) {
  if (val) return std::visit([](auto&& out) -> T& {return *out;}, val.value());
  return wsv;
}

template <class T, class WorkspaceVariable, class VariantT>
T& select_inout(const WorkspaceVariable wsv, std::optional<VariantT>& val) {
  if (val) return std::visit([](auto&& out) -> T& {return *out;}, val.value());
  return wsv;
}

template <class T, class WorkspaceVariable, class VariantT>
const T& select_in(const WorkspaceVariable wsv, const std::optional<VariantT>& val) {
  if (val) return std::visit([](auto&& out) -> const T& {return *out;}, val.value());
  return wsv;
}

template <class SelectT, class VariantT>
SelectT select_wvv(VariantT val) {
  return std::visit([](auto&& out) -> SelectT {return *out;}, val);
}
}  // namespace Python

#endif  // python_interface_h
