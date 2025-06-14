#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <parameters.h>
#include <workspace.h>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "debug.h"
#include "hpy_arts.h"
#include "hpy_vector.h"
#include "nanobind/nanobind.h"
#include "python_interface.h"

extern Parameters parameters;

namespace Python {
Index create_workspace_gin_default_internal(Workspace& ws, const String& key);

std::filesystem::path correct_include_path(
    const std::filesystem::path& path_copy) {
  std::filesystem::path path = path_copy;
  for (auto& prefix : parameters.includepath) {
    if (std::filesystem::is_regular_file(
            path = std::filesystem::path(prefix.c_str()) / path_copy))
      break;
  }

  ARTS_USER_ERROR_IF(not std::filesystem::is_regular_file(path),
                     "Cannot find file: {}\nSearch path(s): {:B,}",
                     path_copy.string(),
                     parameters.includepath)

  // Must add the direcory to include paths as controlfiles know where they are
  std::filesystem::path dir_path = path;
  dir_path.remove_filename();
  String new_inc = dir_path.string();
  if (std::find(parameters.includepath.begin(),
                parameters.includepath.end(),
                new_inc) == parameters.includepath.end())
    parameters.includepath.emplace_back(new_inc);

  return path;
}

void py_agenda(py::module_& m) try {
  py::class_<CallbackOperator> cbd(m, "CallbackOperator");
  workspace_group_interface(cbd);
  cbd.def(
         "__init__",
         [](CallbackOperator* cb,
            const std::function<void(Workspace&)>& f,
            const std::vector<std::string>& i,
            const std::vector<std::string>& o) {
           new (cb) CallbackOperator(
               [f](Workspace& ws) {
                 py::gil_scoped_acquire gil{};
                 f(ws);
               },
               i,
               o);
         },
         "f"_a,
         "inputs"_a  = std::vector<std::string>{},
         "outputs"_a = std::vector<std::string>{},
         "Initialize as structured call")
      .def("__call__", [](CallbackOperator& f, Workspace& ws) { f(ws); });

  py::class_<Wsv>(m, "Wsv")
      .def_prop_ro(
          "value",
          [](Wsv& v) { return to_py(v); },
          "A workspace variable")
      .def(
          "__str__",
          [](py::object& x) { return x.attr("value").attr("__str__")(); },
          py::is_operator())
      .def(
          "__repr__",
          [](py::object& x) { return x.attr("value").attr("__repr__")(); },
          py::is_operator())
      .def(
          "__format__",
          [](py::object& x, py::object& fmt) {
            return x.attr("value").attr("__format__")(fmt);
          },
          py::is_operator())
      .doc() = "A workspace variable wrapper - no manual use required";

  py::class_<Method> methods(m, "Method");
  str_interface(methods);
  methods
      .def(
          "__init__",
          [](Method* me,
             const std::string& n,
             const std::vector<std::string>& a,
             const std::unordered_map<std::string, std::string>& kw) {
            new (me) Method{n, a, kw};
          },
          "name"_a,
          "args"_a,
          "kwargs"_a,
          "A named method with args and kwargs")
      .def(
          "__init__",
          [](Method* me, const std::string& n, const py::object * const v) {
            new (me) Method{n, from(v).copied()};
          },
          "name"_a,
          "wsv"_a,
          "A method that sets a workspace variable")
      .def_prop_ro(
          "val",
          [](const Method& method) -> py::object {
            const auto& x = method.get_setval();
            if (x) return to_py(x.value());
            return py::none();
          },
          "The value (if any) of a set method")
      .def_prop_ro(
          "name",
          [](const Method& method) { return method.get_name(); },
          "The name of the method")
      .doc() = "The method class of ARTS";

  py::class_<Agenda> ag(m, "Agenda");
  workspace_group_interface(ag);
  ag.def(py::init<std::string>(), "name"_a, "Create with name")
      .def("document",
           &Agenda::sphinx_list,
           "prep"_a = std::string_view{"- "},
           "Returns a list of methods")
      .def("add",
           &Agenda::add,
           "method"_a.none(false),
           R"--(
Adds a method to the Agenda

All workspace variables are defaulted, and all GIN with defaults
create anonymous workspace variables.  All input that are not 
workspace variables are added to the workspace

The input order takes priority over the named argument order,
so Copy(a, out=b) will not even see the b variable.
)--")
      .def(
          "execute",
          [](Agenda& a, Workspace& ws) { a.execute(ws); },
          "ws"_a,
          "Executes the agenda on the provided workspace")
      .def(
          "finalize",
          [](Agenda& a, bool fix) { a.finalize(fix); },
          "fix"_a = false,
          "Finalize the agenda, making it possible to use it in the workspace")
      .def_prop_ro(
          "name",
          [](const Agenda& agenda) { return agenda.get_name(); },
          "The name of the agenda")
      .def_prop_ro(
          "methods",
          [](const Agenda& agenda) { return agenda.get_methods(); },
          "The methods of the agenda");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize agendas\n{}", e.what()));
}
}  // namespace Python
