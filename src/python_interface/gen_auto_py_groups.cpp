#include <workspace.h>

#include <iostream>

void groups(const std::string& fname) {
  const auto& wsgs = internal_workspace_groups();

  std::ofstream hos(fname + ".h");

  hos << R"--(#pragma once

#include <auto_wsg.h>
#include <python_interface_value_type.h>

#include <nanobind/nanobind.h>

#include <hpy_vector.h>

)--";

  for (auto& [group, wsg] : wsgs) {
    hos << "NB_MAKE_OPAQUE(Array<" << group << ">);\n";
  }

  std::string_view newline = "\n";
  hos << "namespace Python {\nusing PyWSV = std::variant<";
  for (auto& [group, wsg] : wsgs) {
    hos << std::exchange(newline, ",\n");
    if (wsg.value_type)
      hos << "  ValueHolder<" << group << ">";
    else
      hos << "  std::shared_ptr<" << group << ">";
  }
  hos << ">;\n}\n";

  hos << R"--(
namespace Python {
namespace py = nanobind;

Wsv from(py::object* const x);
Wsv from(const py::object* const x);

std::string type(const py::object * const x);

template <WorkspaceGroup T>
std::string type(const T* const) {
  return std::string{WorkspaceGroupInfo<T>::name};
}

template <WorkspaceGroup T>
std::string type(const ValueHolder<T>* const) {
  return std::string{WorkspaceGroupInfo<T>::name};
}
}  // namespace Python

)--";

  std::ofstream cos(fname + ".cpp");

  cos << R"--(#include <py_auto_wsg.h>

namespace Python {
)--";

  cos << R"--(

Wsv from(const py::object * const x) {
  if (not x or x -> is_none()) throw std::runtime_error("Cannot convert None to workspace variable.");

  return py::cast<Wsv>(*x);
}

Wsv from(py::object * const x) {
  if (not x or x -> is_none()) throw std::runtime_error("Cannot convert None to workspace variable.");

  return py::cast<Wsv>(*x);
}

std::string type(const py::object * const x) {
  if (not x or x -> is_none()) return "NoneType";

  return py::cast<std::string>(py::str(py::type_name(*x)));
}
}  // namespace Python
)--";
}

int main() { groups("py_auto_wsg"); }
