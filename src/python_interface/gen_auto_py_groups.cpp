#include <workspace.h>

#include <iostream>
#include <string>

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
    if (wsg.map_type) {
      hos << "NB_MAKE_OPAQUE(" << group << ");\n";
    }
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
  return std::visit([](auto v) { return from(std::move(v)); }, wsv.value());
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
        unwrap_stars(std::format("{}\n{}{}",
                                 wsg.desc,
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

void agenda_operators() {
  const auto& wsv = internal_workspace_variables();

  std::ofstream cpp("py_auto_agenda_operators.cpp");

  cpp << R"(#include <python_interface.h>

#include <hpy_arts.h>
#include <nanobind/stl/function.h>

namespace Python {

void py_auto_agenda_operators(py::module_& m) {
)";

  for (auto& [name, ag] : internal_workspace_agendas()) {
    std::vector<std::string> input;
    std::vector<std::string> vars;
    std::string params;
    std::string retval;

    input.reserve(ag.input.size());
    vars.reserve(ag.input.size());
    params.reserve(ag.input.size());
    for (auto&& v : ag.input) {
      input.push_back(std::format("const {}& {}", wsv.at(v).type, v));
      vars.push_back(std::format(R"("{}"_a)", v));
      params += std::format(
          R"({0} : :class:`~pyarts.arts.{1}`
     {2} See also :attr:`~pyarts.workspace.Workspace.{0}`.
)",
          v,
          wsv.at(v).type,
          unwrap_stars(short_doc(v)));
    }

    retval.reserve(ag.output.size());
    for (auto&& v : ag.output) {
      retval += std::format(
          R"({0} : :class:`~pyarts.arts.{1}`
     {2} See also :attr:`~pyarts.workspace.Workspace.{0}`.
)",
          v,
          wsv.at(v).type,
          unwrap_stars(short_doc(v)));
    }

    std::print(cpp,
               R"(
  py::class_<{0}Operator> {0}_operator(m, "{0}Operator");
    {0}_operator.def("__init__",
            []({0}Operator* op, {0}Operator::func_t f) {{
              new (op) {0}Operator([f = std::move(f)]({1:,}) {{
                py::gil_scoped_acquire gil{{}};
                return f({2:,});
              }});
            }})
        .def(
            "__call__",
            []({0}Operator& f, {1:,}) {{
              return f.f({2:,});
            }},
            {3:,},
            R"-x-(Execute the method directly in python

Parameters
----------
{4}

Returns
-------
{5}
)-x-"
            );
    workspace_group_interface({0}_operator);
    py::implicitly_convertible<{0}Operator::func_t,
                               {0}Operator>();
    {0}_operator.doc() = R"-x-(This is the operator for free customization of the agenda: :attr:`~pyarts.workspace.Workspace.{0}`.

The python meta-code to make use of this operator class instead of the standard agenda in a workspace reads like this:

.. code-block:: python

    import pyarts

    ws = pyarts.Workspace()

    def my_{0}_operator({2:,}):
        ...
        # custom code that creates or modifies {6:,}
        ...
        return {6:,}
    
    ws.{0}SetOperator(f=my_{0}_operator)

You are free to put whatever custom code you want in the operator.
The types returned must be convertible to the types of the workspace variables being returned.
The input to the operator from the workspace will be the types of the workspace variables
that are passed to the agenda.
Failure to follow these rules will result in a runtime error.

.. note::

      There might be constraints on the input and output passed to the operator.
      Check the documentation of the agenda for more details.

.. warning::

      Accessing global variables in the operator is not thread safe.
      It is very likely ARTS will access your operator in parallel.
      You need to make sure that the operator is thread safe.
)-x-";
)",
               name,
               input,
               ag.input,
               vars,
               params,
               retval,
               ag.output);
  }

  cpp << "}}\n";
}

int main() {
  groups("py_auto_wsg");
  groupdocs("py_auto_wsgdocs.h");
  agenda_operators();
}
