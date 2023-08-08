#ifndef python_interface_h
#define python_interface_h

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

template <WorkspaceGroup T>
T& select_gin(const std::optional<T*>& a, const std::string& name) {
  return a ? **a
           : throw std::runtime_error(
                 var_string("Output variable ", std::quoted(name), " not set"));
}

template <WorkspaceGroup T>
const T& select_gin(const std::optional<T*>& a, const T& defval) {
  return a ? **a : defval;
}

template <typename... Ts>
std::variant<Wsv, py::object *> select_wsv(
    const std::optional<std::variant<Ts*...>>& a, const std::string& name) {
  return a ? std::visit(
                 [](auto& x) {
                   if constexpr (WorkspaceGroup<
                                     std::remove_cvref_t<decltype(x)>>)
                     return Wsv(x);
                   else
                     return x;
                 },
                 *a)
           : throw std::runtime_error(
                 var_string("Variable ", std::quoted(name), " not set"));
}
}  // namespace Python

#endif  // python_interface_h
