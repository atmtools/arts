#ifndef python_interface_h
#define python_interface_h

#include <py_auto_wsg.h>
#include <pybind11/attr.h>
#include <pybind11/cast.h>
#include <pybind11/chrono.h>
#include <pybind11/complex.h>
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

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include "auto_wsg.h"

//! Contains a bunch of helper functions to manipulate python objects inside C++
namespace Python {
namespace py = pybind11;

template <typename T, typename ... Ts>
using artsclass = py::class_<T, Ts..., std::shared_ptr<T>>;

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
T& select_out(std::optional<std::shared_ptr<U>>& x, Workspace& ws, const std::string& name) {
  if (x) return *x.value();
  return ws.get_or<T>(name);
}

template <WorkspaceGroup T>
T& select_out(std::optional<ValueHolder<T>>& x, Workspace& ws, const std::string& name) {
  if (x) return *x;
  return ws.get_or<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
T& select_gout(std::optional<std::shared_ptr<U>>& x, const std::string& name) {
  if (x) return *x.value();
  throw std::runtime_error(var_string("Unknown ouput: ", std::quoted(name)));
}

template <WorkspaceGroup T>
T& select_gout(std::optional<ValueHolder<T>>& x, const std::string& name) {
  if (x) return x.value();
  throw std::runtime_error(var_string("Unknown ouput: ", std::quoted(name)));
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
T& select_inout(std::optional<std::shared_ptr<U>>& x, const Workspace& ws, const std::string& name) {
  if (x) return *x.value();
  return ws.get<T>(name);
}

template <WorkspaceGroup T>
T& select_inout(std::optional<ValueHolder<T>>& x, const Workspace& ws, const std::string& name) {
  if (x) return *x;
  return ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
T& select_in(const std::optional<std::shared_ptr<U>>& x, const Workspace& ws, const std::string& name) {
  if (x) return *x.value();
  return ws.get<T>(name);
}

template <WorkspaceGroup T>
const T& select_in(const std::optional<ValueHolder<T>>& x, const Workspace& ws, const std::string& name) {
  if (x) return x.value();
  return ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
T& select_gin(const std::optional<std::shared_ptr<U>>& x, const std::string& name) {
  if (x) return *x.value();
  throw std::runtime_error(var_string("Unknown input: ", std::quoted(name)));
}

template <WorkspaceGroup T>
const T& select_gin(const std::optional<ValueHolder<T>>& x, const std::string& name) {
  if (x) return x.value();
  throw std::runtime_error(var_string("Unknown input: ", std::quoted(name)));
}

template <WorkspaceGroup T, PythonWorkspaceGroup U>
T& select_gin(const std::optional<std::shared_ptr<U>>& x, T& defval) {
  if (x) return *x.value();
  return defval;
}

template <WorkspaceGroup T>
const T& select_gin(const std::optional<ValueHolder<T>>& x, const T& defval) {
  if (x) return x.value();
  return defval;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
}  // namespace Python

#endif  // python_interface_h
