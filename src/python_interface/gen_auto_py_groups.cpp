#include <workspace.h>
#include <workspace_group_friends.h>

#include <iostream>
#include <ranges>
#include <string>

#include "pydocs.h"

void implement_convert_const_py_object() {
  const auto& wsgs = internal_workspace_groups();

  std::ofstream os("py_auto_wsg_convert.cpp");
  std::println(os, R"--(#include <py_auto_wsg.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/nanobind.h>

#include "py_auto_options.h"

namespace Python {{
bool convert_ref(Wsv& wsv, const py::object * const x) {{
  py::gil_scoped_acquire gil{{}};
  if (not x or x -> is_none()) throw std::runtime_error("Cannot convert None to workspace variable.");
)--");

  for (auto& [group, wsg] : wsgs) {
    if (wsg.value_type) {
      if (group == "Numeric"sv) {
        std::println(os,
                     R"(  if (py::isinstance<ValueHolder<Numeric>>(*x)) {{
    wsv = Numeric{{py::float_(*x)}};
    return true;
  }})");
      } else if (group == "String"sv) {
        std::println(os,
                     R"( if (py::isinstance<ValueHolder<String>>(*x)) {{
    wsv = String{{py::str(*x).c_str()}};
    return true;
  }})");
      } else if (group == "Index") {
        std::println(os,
                     R"(  if (py::isinstance<ValueHolder<Index>>(*x)) {{
    wsv = Index{{py::int_(*x)}};
    return true;
  }})");
      } else {
        std::print(os,
                   "static_assert(false, \"Missing custom code for {0}\");",
                   group);
      }
    } else {
      std::println(os,
                   R"(  if (py::isinstance<{0}>(*x)) {{
    wsv = py::cast<std::shared_ptr<{0}>>(*x, false);
    return true;
  }})",
                   group);
    }
  }

  os << "  return false;}\n\n";
  std::println(os,
               R"--(bool convert_cast(Wsv& wsv, const py::object * const x) {{
  py::gil_scoped_acquire gil{{}};
  if (not x or x -> is_none()) throw std::runtime_error("Cannot convert None to workspace variable.");
)--");

  for (auto& [group, wsg] : wsgs) {
    if (wsg.value_type) {
      std::print(os,
                 R"(
  ValueHolder<{0}> _val{0}{{}};
  if (py::try_cast(*x, _val{0}, true)) {{
    wsv = std::move(_val{0}.val);
    return true;
  }})",
                 group);
    } else {
      std::print(os,
                 R"(
  {0} _val{0}{{}};
  if (py::try_cast(*x, _val{0}, true)) {{
    wsv = _val{0};
    return true;
  }})",
                 group);
    }
  }

  os << "\n  return false;\n}\n}\n";
}

void implement_from_const_py_object() {
  std::ofstream os("py_auto_wsg_from_const_py_object.cpp");
  std::println(os, R"--(#include <py_auto_wsg.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/nanobind.h>

#include "py_auto_options.h"

namespace Python {{
Wsv from(const py::object * const x) {{
  Wsv wsv{{}};

  if (convert_ref(wsv, x)) return wsv;
  if (convert_cast(wsv, x)) return wsv;

  throw std::runtime_error("Cannot convert pure python object to workspace variable.");
}}
}}  // namespace Python
)--");
}

void implement_from_py_object() {
  const auto& wsgs = internal_workspace_groups();

  std::ofstream os("py_auto_wsg_from_py_object.cpp");
  os << R"--(#include <py_auto_wsg.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/nanobind.h>

#include "py_auto_options.h"

namespace Python {
)--";
  os << R"--(

Wsv from(py::object * const x) {
  py::gil_scoped_acquire gil{};
  if (not x or x -> is_none()) throw std::runtime_error("Cannot have None as workspace variable.");
)--";

  for (auto& [group, wsg] : wsgs) {
    if (wsg.value_type) {
      os << "  if (py::isinstance<ValueHolder<" << group
         << ">>(*x)) return py::cast<ValueHolder<" << group
         << ">>(py::object(x->attr(\"value\")), false).val;\n";
    } else {
      os << "  if (py::isinstance<" << group
         << ">(*x)) return py::cast<std::shared_ptr<" << group
         << ">>(*x, false);\n";
    }
  }

  os << R"--(
  
  throw std::runtime_error("Cannot convert pure python object to workspace variable.");
}
}  // namespace Python
)--";
}

void implement_string_type() {
  std::ofstream os("py_auto_wsg_string_type.cpp");
  os << R"--(#include <py_auto_wsg.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/nanobind.h>

#include "py_auto_options.h"

namespace Python {
std::string type(const py::object * const x) {
  py::gil_scoped_acquire gil{};
  if (not x or x -> is_none()) return "NoneType";

  return py::cast<std::string>(py::str(py::type_name(*x)));
}
}  // namespace Python
)--";
}

void implement_to_py_wsv() {
  const auto& wsgs = internal_workspace_groups();

  std::ofstream os("py_auto_wsg_to_py_wsv.cpp");

  std::println(os, R"--(#include <py_auto_wsg.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/nanobind.h>

#include "py_auto_options.h"

namespace Python {{
py::object to_py(const Wsv& wsv) {{
  switch (wsv.value_index()) {{)--");

  for (auto& [group, wsg] : wsgs) {
    if (wsg.value_type) {
      std::println(
          os,
          "    case WorkspaceGroupInfo<{0}>::index:\n      return py::cast<ValueHolder<{0}>>(wsv.share<{0}>());",
          group);
    } else {
      std::println(
          os,
          "    case WorkspaceGroupInfo<{0}>::index:\n      return py::cast<std::shared_ptr<{0}>>(wsv.share<{0}>());",
          group);
    }
  }

  os << R"--(  }
  return py::none();
}
}  // namespace Python
)--";
}

void groups(const std::string& fname) {
  const auto& wsgs = internal_workspace_groups();

  implement_from_const_py_object();
  implement_from_py_object();
  implement_string_type();
  implement_convert_const_py_object();
  implement_to_py_wsv();

  std::ofstream hos(fname + ".h");

  hos << R"--(#pragma once

#include <auto_wsg.h>
#include <python_interface_value_type.h>

#include <nanobind/nanobind.h>

#include <hpy_vector.h>

)--";

  for (auto& [group, wsg] :
       std::array{wsgs, workspace_group_friends()} | stdv::join) {
    hos << "NB_MAKE_OPAQUE(Array<" << group << ">);\n";
    if (wsg.map_type) hos << "NB_MAKE_OPAQUE(" << group << ");\n";
  }

  hos << R"--(
namespace Python {
namespace py = nanobind;

bool convert_ref(Wsv& wsv, const py::object * const x);
bool convert_cast(Wsv& wsv, const py::object * const x);
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
Wsv from_py(std::shared_ptr<T> wsv) {
  return wsv;
}

template <WorkspaceGroup T>
Wsv from_py(ValueHolder<T> wsv) {
  std::shared_ptr<T> copy = wsv.val;
  return copy;
}

py::object to_py(const Wsv& wsv);

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
}

void groupdocs(const std::string& fname) {
  const auto& wsgs = internal_workspace_groups();

  std::ofstream osh(fname + ".h");
  std::ofstream osc(fname + ".cpp");

  std::print(osh, R"(#pragma once

#include <workspace.h>

template <typename T>
struct PythonWorkspaceGroupInfo {{static std::string_view desc() = delete;}};
)");

  std::println(osc, R"(#include "py_auto_wsgdocs.h"
)");

  for (auto& [group, wsg] :
       std::array{wsgs, workspace_group_friends()} | stdv::join) {
    const auto info =
        unwrap_stars(std::format("{}\n{}{}",
                                 wsg.desc,
                                 Python::group_generics_inout(group),
                                 Python::group_workspace_types(group)));

    std::println(osh,

                 R"(template <>
struct PythonWorkspaceGroupInfo<{0}> {{
  static std::string_view desc();
}};
)",
                 group);

    std::println(osc,
                 R"(std::string_view PythonWorkspaceGroupInfo<{0}>::desc() {{
  return
  R"-PYARTSDOCSAUTO-({1})-PYARTSDOCSAUTO-"sv;
}}
)",
                 group,
                 info);
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
          R"({0} : :class:`~pyarts3.arts.{1}`
     {2} See also :attr:`~pyarts3.workspace.Workspace.{0}`.
)",
          v,
          wsv.at(v).type,
          unwrap_stars(short_doc(v)));
    }

    retval.reserve(ag.output.size());
    for (auto&& v : ag.output) {
      retval += std::format(
          R"({0} : :class:`~pyarts3.arts.{1}`
     {2} See also :attr:`~pyarts3.workspace.Workspace.{0}`.
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
    generic_interface({0}_operator);
    py::implicitly_convertible<{0}Operator::func_t,
                               {0}Operator>();
    {0}_operator.doc() = R"-x-(This is the operator for free customization of the agenda: :attr:`~pyarts3.workspace.Workspace.{0}`.

The python meta-code to make use of this operator class instead of the standard agenda in a workspace reads like this:

.. code-block:: python

    import pyarts3 as pyarts

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
  groupdocs("py_auto_wsgdocs");
  agenda_operators();
}
