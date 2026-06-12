#ifndef python_interface_h
#define python_interface_h

#include <nanobind/nanobind.h>
#include <py_auto_options.h>
#include <py_auto_wsg.h>
#include <workspace.h>

#include <memory>
#include <variant>

#include "hpy_opaque.h"

using ssize_t = Py_ssize_t;

//! Contains a bunch of helper functions to manipulate python objects inside C++
namespace Python {
namespace py = nanobind;
using namespace py::literals;

template <typename... T>
bool has_selected_value(const std::variant<std::shared_ptr<T>...>& x) {
  return std::visit([](const auto& ptr) { return static_cast<bool>(ptr); }, x);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
T& select_out(T* const x, Workspace& ws, const char* const name) {
  return x ? *x : ws.get_or<T>(name);
}

template <WorkspaceGroup T>
T& select_out(ValueHolder<T>* const x, Workspace& ws, const char* const name) {
  return x ? static_cast<T&>(*x) : ws.get_or<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
T& select_gout(T* const x, Workspace& ws, const char* const name) {
  return x ? *x : ws.get_or<T>(name);
}

template <WorkspaceGroup T>
T& select_gout(ValueHolder<T>* const x, Workspace& ws, const char* const name) {
  return x ? static_cast<T&>(*x) : ws.get_or<T>(name);
}

template <typename... T>
std::variant<std::shared_ptr<T>...>& select_gout(
    std::variant<std::shared_ptr<T>...>* const x,
    Workspace&,
    const char* const name) {
  if (not x or not has_selected_value(*x)) {
    throw std::runtime_error(std::format("Unknown output: \"{}\"", name));
  }

  return *x;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
T& select_inout(T* const x, const Workspace& ws, const char* const name) {
  return x ? *x : ws.get<T>(name);
}

template <WorkspaceGroup T>
T& select_inout(ValueHolder<T>* const x,
                const Workspace& ws,
                const char* const name) {
  return x ? static_cast<T&>(*x) : ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
const T& select_in(const T* const x,
                   const Workspace& ws,
                   const char* const name) {
  return x ? *x : ws.get<T>(name);
}

template <WorkspaceGroup T>
const T& select_in(const ValueHolder<T>* const x,
                   const Workspace& ws,
                   const char* const name) {
  return x ? static_cast<const T&>(*x) : ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
const T& select_gin(const T* const x, const char* const name) {
  return x ? *x
           : throw std::runtime_error(
                 std::format("Unknown input: \"{}\"", name));
}

template <WorkspaceGroup T>
const T& select_gin(const ValueHolder<T>* const x, const char* const name) {
  return x ? static_cast<const T&>(*x)
           : throw std::runtime_error(
                 std::format("Unknown input: \"{}\"", name));
}

template <typename... T>
const std::variant<std::shared_ptr<T>...>& select_gin(
    const std::variant<std::shared_ptr<T>...>* const x,
    const char* const name) {
  if (not x or not has_selected_value(*x)) {
    throw std::runtime_error(std::format("Unknown input: \"{}\"", name));
  }

  return *x;
}

template <WorkspaceGroup T>
const T& select_gin(const T* const x, const T& defval) {
  return x ? *x : defval;
}

template <WorkspaceGroup T>
const T& select_gin(const ValueHolder<T>* const x, const T& defval) {
  return x ? static_cast<const T&>(*x) : defval;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
}  // namespace Python

template <typename T>
struct std::hash<Python::ValueHolder<T>> {
  static std::size_t operator()(const Python::ValueHolder<T>& x) {
    return std::hash<T>{}(x->val);
  }
};

#endif  // python_interface_h
