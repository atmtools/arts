#ifndef python_interface_h
#define python_interface_h

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
#include <variant>

#include "callback.h"

//! Contains a bunch of helper functions to manipulate python objects inside C++
namespace Python {
namespace py = pybind11;

//! To suppress warnings about unused variables
void print(WorkspaceGroup auto& x) {py::print(x);}

template <WorkspaceGroup T>
T& select_out(std::optional<T*>& a, Workspace& ws, const std::string& name) {
  return a ? **a : ws.get_or<T>(name);
}

template <WorkspaceGroup T>
T& select_inout(std::optional<T*>& a,
                const Workspace& ws,
                const std::string& name) {
  return a ? **a : ws.get<T>(name);
}

template <WorkspaceGroup T>
const T& select_in(const std::optional<T*>& a,
                   const Workspace& ws,
                   const std::string& name) {
  return a ? **a : ws.get<T>(name);
}

template <WorkspaceGroup T>
T& select_gout(std::optional<T*>& a, const std::string& name) {
  return a ? **a
           : throw std::runtime_error(
                 var_string("Output variable ", std::quoted(name), " not set"));
}

template <WorkspaceGroup ... Ts>
std::variant<Ts*...>& select_gout(std::optional<std::variant<Ts*...>>& a, const std::string& name) {
  return a ? *a
           : throw std::runtime_error(
                 var_string("Output variable ", std::quoted(name), " not set"));
}

inline Wsv& select_gout(std::optional<Wsv *> & a, const std::string& name) {
  return a ? **a
           : throw std::runtime_error(
                 var_string("Output variable ", std::quoted(name), " not set"));
}

template <WorkspaceGroup T>
const T& select_gin(const std::optional<T*>& a, const std::string& name) {
  return a ? **a
           : throw std::runtime_error(
                 var_string("Input variable ", std::quoted(name), " not set"));
}

template <WorkspaceGroup ... Ts>
const std::variant<Ts*...>& select_gin(const std::optional<std::variant<Ts*...>>& a, const std::string& name) {
  return a ? *a
           : throw std::runtime_error(
                 var_string("Input variable ", std::quoted(name), " not set"));
}

inline const Wsv& select_gin(const std::optional<Wsv *> & a, const std::string& name) {
  return a ? **a
           : throw std::runtime_error(
                 var_string("Input variable ", std::quoted(name), " not set"));
}

template <WorkspaceGroup T>
const T& select_gin(const std::optional<T*>& a, const T& defval) {
  return a ? **a : defval;
}
}  // namespace Python

#endif  // python_interface_h
