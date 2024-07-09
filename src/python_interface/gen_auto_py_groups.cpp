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

template <WorkspaceGroup T>
PyWSV from(std::shared_ptr<T> wsv) {
  if constexpr (WorkspaceGroupInfo<T>::value_type)
    return ValueHolder<T>(std::move(wsv));
  else
    return wsv;
}

PyWSV from(const std::shared_ptr<Wsv>& wsv);

template <WorkspaceGroup T>
Wsv from_py(std::shared_ptr<T> wsv) {
  return wsv;
}

template <WorkspaceGroup T>
Wsv from_py(ValueHolder<T> wsv) {
  std::shared_ptr<T> copy = wsv.val;
  return copy;
}

Wsv from_py(const PyWSV& wsv);
}  // namespace Python

)--";

  std::ofstream cos(fname + ".cpp");

  cos << R"--(#include <py_auto_wsg.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/nanobind.h>

NB_MAKE_OPAQUE(SpeciesEnum)

namespace Python {
)--";

  cos << R"--(

Wsv from(const py::object * const x) {
  if (not x or x -> is_none()) throw std::runtime_error("Cannot convert None to workspace variable.");
  
)--";

  for (auto& [group, wsg] : wsgs) {
    if (wsg.value_type) {
      cos << "  if (py::isinstance<ValueHolder<" << group
          << ">>(*x)) return py::cast<ValueHolder<" << group
          << ">>(py::object(x->attr(\"value\")), false).val;\n";
    } else {
      cos << "  if (py::isinstance<" << group
          << ">(*x)) return py::cast<std::shared_ptr<" << group
          << ">>(*x, false);\n";
    }
  }

  cos << R"--(
  return from_py(py::cast<PyWSV>(*x));
}

Wsv from(py::object * const x) {
  if (not x or x -> is_none()) throw std::runtime_error("Cannot have None as workspace variable.");
)--";

  for (auto& [group, wsg] : wsgs) {
    if (wsg.value_type) {
      cos << "  if (py::isinstance<ValueHolder<" << group
          << ">>(*x)) return py::cast<ValueHolder<" << group
          << ">>(py::object(x->attr(\"value\")), false).val;\n";
    } else {
      cos << "  if (py::isinstance<" << group
          << ">(*x)) return py::cast<std::shared_ptr<" << group
          << ">>(*x, false);\n";
    }
  }

  cos << R"--(
  
  throw std::runtime_error("Cannot convert pure python object to workspace variable.");
}

std::string type(const py::object * const x) {
  if (not x or x -> is_none()) return "NoneType";

  return py::cast<std::string>(py::str(py::type_name(*x)));
}

PyWSV from(const std::shared_ptr<Wsv>& wsv) {
  return std::visit([](auto v) { return from(std::move(v)); }, wsv->value);
}

Wsv from_py(const PyWSV& wsv) {
  return std::visit([](auto v) { return from_py(std::move(v)); }, wsv);
}
}  // namespace Python
)--";
}

int main() { groups("py_auto_wsg"); }
