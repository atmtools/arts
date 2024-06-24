#ifndef python_interface_h
#define python_interface_h

#include <nanobind/nanobind.h>
#include <py_auto_wsg.h>
#include <workspace.h>

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <ios>
#include <istream>
#include <iterator>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <streambuf>
#include <type_traits>
#include <variant>

#include "auto_wsg.h"
#include "enums.h"

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
T& select_gout(T* const x, const char* const name) {
  return x ? *x
           : throw std::runtime_error(
                 var_string("Unknown ouput: ", std::quoted(name)));
}

template <WorkspaceGroup T>
T& select_gout(ValueHolder<T>* const x, const char* const name) {
  return x ? static_cast<T&>(*x)
           : throw std::runtime_error(
                 var_string("Unknown ouput: ", std::quoted(name)));
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
                 var_string("Unknown input: ", std::quoted(name)));
}

template <WorkspaceGroup T>
const T& select_gin(const ValueHolder<T>* const x, const char* const name) {
  return x ? static_cast<const T&>(*x)
           : throw std::runtime_error(
                 var_string("Unknown input: ", std::quoted(name)));
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

#endif  // python_interface_h
