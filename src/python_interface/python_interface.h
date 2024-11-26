#ifndef python_interface_h
#define python_interface_h

#include <nanobind/nanobind.h>
#include <py_auto_options.h>
#include <py_auto_wsg.h>
#include <workspace.h>

#include "workspace_class.h"

using ssize_t = Py_ssize_t;

//! Contains a bunch of helper functions to manipulate python objects inside C++
namespace Python {
namespace py = nanobind;
using namespace py::literals;

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
  std::size_t operator()(const Python::ValueHolder<T>& x) const {
    return std::hash<T>{}(x->val);
  }
};

#endif  // python_interface_h
