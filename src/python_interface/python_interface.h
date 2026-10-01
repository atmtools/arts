#ifndef python_interface_h
#define python_interface_h

#include <nanobind/nanobind.h>
#include <py_auto_options.h>
#include <py_auto_wsg.h>

// NOTE: this header is included by nearly every file in python_interface/,
// so it deliberately does NOT pull in the full <workspace.h> (which drags in
// the ~9700-line generated auto_wsm.h method declarations and the ~1000-line
// auto_wsa.h agenda declarations neither of which anything below needs).
// Only the Workspace class itself and the WorkspaceGroup type universe
// (already available via py_auto_wsg.h -> auto_wsg.h) are required here.
// Files that actually bind workspace methods/agendas (py_workspace.cpp,
// py_agenda.cpp, py_module.cpp, and the generated py_auto_wsm_*.cpp/
// py_auto_wsa.cpp shards) include <workspace.h> themselves.
#include <workspace_class.h>

#include <memory>
#include <variant>

#include "hpy_opaque.h"

using ssize_t = Py_ssize_t;

//! Contains a bunch of helper functions to manipulate python objects inside C++
namespace Python {
namespace py = nanobind;
using namespace py::literals;

template <WorkspaceGroup T> T& select_out(T* const x, Workspace& ws, const char* const name) {
  return x ? *x : ws.get_or<T>(name);
}

template <WorkspaceGroup T> T& select_out(ValueHolder<T>* const x, Workspace& ws, const char* const name) {
  return x ? static_cast<T&>(*x) : ws.get_or<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T> T& select_gout(T* const x, Workspace& ws, const char* const name) {
  return x ? *x : ws.get_or<T>(name);
}

template <WorkspaceGroup T> T& select_gout(ValueHolder<T>* const x, Workspace& ws, const char* const name) {
  return x ? static_cast<T&>(*x) : ws.get_or<T>(name);
}

template <WorkspaceGroup T> T& select_inout(T* const x, const Workspace& ws, const char* const name) {
  return x ? *x : ws.get<T>(name);
}

template <WorkspaceGroup T> T& select_inout(ValueHolder<T>* const x, const Workspace& ws, const char* const name) {
  return x ? static_cast<T&>(*x) : ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T> const T& select_in(const T* const x, const Workspace& ws, const char* const name) {
  return x ? *x : ws.get<T>(name);
}

template <WorkspaceGroup T>
const T& select_in(const ValueHolder<T>* const x, const Workspace& ws, const char* const name) {
  return x ? static_cast<const T&>(*x) : ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T> const T& select_gin(const T* const x, const char* const name) {
  return x ? *x : throw std::runtime_error(std::format("Unknown input: \"{}\"", name));
}

template <WorkspaceGroup T> const T& select_gin(const ValueHolder<T>* const x, const char* const name) {
  return x ? static_cast<const T&>(*x) : throw std::runtime_error(std::format("Unknown input: \"{}\"", name));
}

template <WorkspaceGroup T> const T& select_gin(const T* const x, const T& defval) { return x ? *x : defval; }

template <WorkspaceGroup T> const T& select_gin(const ValueHolder<T>* const x, const T& defval) {
  return x ? static_cast<const T&>(*x) : defval;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
}  // namespace Python

template <typename T> struct std::hash<Python::ValueHolder<T>> {
  static std::size_t operator()(const Python::ValueHolder<T>& x) { return std::hash<T>{}(x->val); }
};

#endif  // python_interface_h
