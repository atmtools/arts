#include <workspace.h>

#include <iostream>

#include "pydocs.h"

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

PyWSV from(const Wsv& wsv);

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

template <typename T>
concept SharedFromPyable = requires(const std::shared_ptr<T>& x) {
  { from_py(x) } -> std::same_as<Wsv>;
} or requires(const std::shared_ptr<T>& x) {
  { from_py(*x) } -> std::same_as<Wsv>;
};

template <SharedFromPyable ... T>
Wsv from(const std::variant<std::shared_ptr<T>...> * const x)  {
  if (not x) throw std::runtime_error("Cannot convert None to workspace variable.");

  return std::visit([]<typename U>(const std::shared_ptr<U>& y) -> Wsv {
    if constexpr (WorkspaceGroup<U>) return from_py(y);
    else return from_py(*y);
  }, *x);
}

template <SharedFromPyable ... T>
std::string type(const std::variant<std::shared_ptr<T>...> * const x)  {
  if (not x) throw std::runtime_error("Cannot convert None to workspace variable.");

  return std::visit([]<typename U>(const std::shared_ptr<U>& y) -> std::string {
    return type(y.get());
  }, *x);
}
}  // namespace Python

)--";

  std::ofstream cos(fname + ".cpp");

  cos << R"--(#include <py_auto_wsg.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/nanobind.h>

#include "py_auto_options.h"

namespace Python {
)--";

  cos << R"--(

Wsv from(const py::object * const x) {
  if (not x or x -> is_none()) throw std::runtime_error("Cannot convert None to workspace variable.");
  py::gil_scoped_acquire gil{};
  
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
  py::gil_scoped_acquire gil{};
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
  py::gil_scoped_acquire gil{};

  return py::cast<std::string>(py::str(py::type_name(*x)));
}

PyWSV from(const Wsv& wsv) {
  return std::visit([](auto v) { return from(std::move(v)); }, wsv.value);
}

Wsv from_py(const PyWSV& wsv) {
  return std::visit([](auto v) { return from_py(std::move(v)); }, wsv);
}
}  // namespace Python
)--";
}

void groupdocs(const std::string& fname) {
  const auto& wsgs = internal_workspace_groups();

  std::ofstream os(fname);

  os << R"(#pragma once

#include <workspace.h>

template <typename T>
struct PythonWorkspaceGroupInfo {constexpr static const char* desc = "Unknown";};
)";

  for (auto& [group, wsg] : wsgs) {
    const auto info =
        unwrap_stars(var_string(wsg.desc,
                                '\n',
                                Python::group_generics_inout(group),
                                Python::group_workspace_types(group)));

    os << R"(
template <>
struct PythonWorkspaceGroupInfo<)"
       << group << R"(> {
  constexpr static const char* desc = R"-X-()"
       << info << R"()-X-";
};
)";
  }
}

int main() {
  groups("py_auto_wsg");
  groupdocs("py_auto_wsgdocs.h");
}
