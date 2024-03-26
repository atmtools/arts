#ifndef python_interface_h
#define python_interface_h

#include <py_auto_wsg.h>
#include <py_auto_wsg_init.h>
#include <pybind11/attr.h>
#include <pybind11/cast.h>
#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/detail/common.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <pybind11/stl_bind.h>
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
#include "python_interface_groups.h"

using ssize_t = Py_ssize_t;

//! Contains a bunch of helper functions to manipulate python objects inside C++
namespace Python {
namespace py = pybind11;

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
T& select_out(std::optional<U* const>& x,
              Workspace& ws,
              const char* const name) {
  return (x and x.value()) ? *x.value() : ws.get_or<T>(name);
}

template <WorkspaceGroup T>
T& select_out(std::optional<ValueHolder<T>* const>& x,
              Workspace& ws,
              const char* const name) {
  return (x and x.value()) ? static_cast<T&>(*x.value()) : ws.get_or<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
T& select_gout(std::optional<U* const>& x, const char* const name) {
  return (x and x.value()) ? *x.value()
                           : throw std::runtime_error(var_string(
                                 "Unknown ouput: ", std::quoted(name)));
}

template <WorkspaceGroup T>
T& select_gout(std::optional<ValueHolder<T>* const>& x,
               const char* const name) {
  return (x and x.value()) ? *x.value()
                           : throw std::runtime_error(var_string(
                                 "Unknown ouput: ", std::quoted(name)));
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
T& select_inout(std::optional<U* const> x,
                const Workspace& ws,
                const char* const name) {
  return (x and x.value()) ? *x.value() : ws.get<T>(name);
}

template <WorkspaceGroup T>
T& select_inout(std::optional<ValueHolder<T>* const>& x,
                const Workspace& ws,
                const char* const name) {
  return (x and x.value()) ? static_cast<T&>(*x.value()) : ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
const T& select_in(const std::optional<const U* const>& x,
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

template <WorkspaceGroup T, PythonWorkspaceGroup U>
const T& select_gin(const std::optional<const U* const>& x,
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

template <WorkspaceGroup T, PythonWorkspaceGroup U>
const T& select_gin(const std::optional<const U* const>& x, const T& defval) {
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
