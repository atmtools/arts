#ifndef python_interface_h
#define python_interface_h

#include <py_auto_wsg.h>
#include <nanobind/nanobind.h>
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
using namespace nanobind::literals;

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
T& select_out(std::optional< T* const>& x,
              Workspace& ws,
              const char* const name) {
  return x ? *x.value() : ws.get_or<T>(name);
}

template <WorkspaceGroup T>
T& select_out(std::optional<ValueHolder<T>* const>& x,
              Workspace& ws,
              const char* const name) {
  return x ? static_cast<T&>(x) : ws.get_or<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
T& select_gout(std::optional<T* const>& x, const char* const name) {
  return x ? *x.value()
           : throw std::runtime_error(
                 var_string("Unknown ouput: ", std::quoted(name)));
}

template <WorkspaceGroup T>
T& select_gout(std::optional<ValueHolder<T>* const>& x,
               const char* const name) {
  return x ? static_cast<T&>(x)
           : throw std::runtime_error(
                 var_string("Unknown ouput: ", std::quoted(name)));
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
T& select_inout(std::optional<T* const> x,
                const Workspace& ws,
                const char* const name) {
  return x ? *x.value() : ws.get<T>(name);
}

template <WorkspaceGroup T>
T& select_inout(std::optional<ValueHolder<T>* const>& x,
                const Workspace& ws,
                const char* const name) {
  return x ? static_cast<T&>(x) : ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
const T& select_in(const std::optional<const T* const>& x,
                   const Workspace& ws,
                   const char* const name) {
  return (x and x.value()) ? *x.value() : ws.get<T>(name);
}

template <WorkspaceGroup T>
const T& select_in(const std::optional<const ValueHolder<T>* const>& x,
                   const Workspace& ws,
                   const char* const name) {
  return (x and x.value()) ? static_cast<const T&>(*x.value())
                           : ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
const T& select_gin(const std::optional<const T* const>& x,
                    const char* const name) {
  return (x and x.value()) ? *x.value()
                           : throw std::runtime_error(var_string(
                                 "Unknown input: ", std::quoted(name)));
}

template <WorkspaceGroup T>
const T& select_gin(const std::optional<const ValueHolder<T>* const>& x,
                    const char* const name) {
  return (x and x.value()) ? *x.value()
                           : throw std::runtime_error(var_string(
                                 "Unknown input: ", std::quoted(name)));
}

template <WorkspaceGroup T>
const T& select_gin(const std::optional<const T* const>& x, const T& defval) {
  return (x and x.value()) ? *x.value() : defval;
}

template <WorkspaceGroup T>
const T& select_gin(const std::optional<const ValueHolder<T>* const>& x,
                    const T& defval) {
  return (x and x.value()) ? static_cast<const T&>(*x.value()) : defval;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
}  // namespace Python

#endif  // python_interface_h
