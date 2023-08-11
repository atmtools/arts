#include <algorithm>
#include <parameters.h>
#include <memory>
#include <type_traits>

#include "debug.h"
#include "py_macros.h"
#include "python_interface.h"
#include "workspace_agenda_class.h"

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
                     "Cannot find file: ",
                     path_copy,
                     '\n',
                     "Search path: ",
                     parameters.includepath)

  // Must add the direcory to include paths as controlfiles know where they are
  std::filesystem::path dir_path = path;
  dir_path.remove_filename();
  String new_inc = dir_path.c_str();
  if (std::find(parameters.includepath.begin(),
                parameters.includepath.end(),
                new_inc) == parameters.includepath.end())
    parameters.includepath.emplace_back(new_inc);

  return path;
}

void py_agenda(py::module_& m) try {
  artsclass<CallbackFunction>(m, "CallbackFunction")
      .def(py::init([]() { return std::make_shared<CallbackFunction>(); }), py::doc("Initialize as empty call"))
      .PythonInterfaceCopyValue(CallbackFunction)
      .def(py::init([](const std::function<void(const Workspace&)>& f) {
        return std::make_shared<CallbackFunction>(f);
      }), py::doc("Initialize from python callable"))
      .def(
          "__call__",
          [](CallbackFunction& f, Workspace& ws) { f(ws); })
      .PythonInterfaceWorkspaceDocumentationExtra(CallbackFunction, R"(

An instance of this object can be called taking only an
instance of the workspace as input.

Example
-------
>>> import pyarts
>>> def print_ws(ws):
>>>     print(ws.iy_unit.value)
>>> ws = pyarts.workspace.Workspace()
>>> cb = pyarts.arts.CallbackFunction(print_ws)
>>> cb(ws)
1
)");
  py::implicitly_convertible<std::function<void(Workspace&)>,
                             CallbackFunction>();

  artsclass<Agenda>(m, "Agenda")
      .def(py::init([]() { return std::make_shared<Agenda>(); }), "Create empty")
      .PythonInterfaceWorkspaceVariableConversion(Agenda)
      .PythonInterfaceFileIO(Agenda)
      .def(
          "add_workspace_method",
          [](Agenda& a,
             const char* name,
             const py::args& args,
             const py::kwargs& kwargs) {
            //pass for now
          },
          py::arg("name").none(false),
          R"--(
Adds a named method to the Agenda

All workspace variables are defaulted, and all GIN with defaults
create anonymous workspace variables.  All input that are not 
workspace variables are added to the workspace

The input order takes priority over the named argument order,
so Copy(a, out=b) will not even see the b variable.
)--")
      .def(
          "add_callback_method",
          [](Agenda& a, const CallbackFunction& f) {
            a.add(Method("@", f, true));
            a.add(Method("CallbackFunctionExecute", {"@"}));
          },
          py::doc(R"--(
Adds a callback method to the Agenda

Parameters
----------
    f : CallbackFunction
        A method that takes ws and only ws as input
)--"),
          py::arg("f"))
      .def(
          "execute",
          [](Agenda& a, Workspace& ws) {
            a.execute(ws);
          },
          "Executes the agenda as if it was the main agenda")
      .def(
          "finalize",
          [](Agenda& a) {
            a.finalize();
          })
      .def("__repr__",
           [](Agenda& a) { return var_string("Agenda ", a.get_name()); })
      .def("__str__",
           [](Agenda& a) {
             return var_string(a);
           })
      .PythonInterfaceWorkspaceDocumentation(Agenda);

  artsclass<ArrayOfAgenda>(m, "ArrayOfAgenda")
      .def(py::init([]() { return std::make_shared<ArrayOfAgenda>(); }), "Create empty")
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfAgenda)
      .PythonInterfaceCopyValue(ArrayOfAgenda)
      .def(py::init([](std::vector<Agenda> va) {
        for (auto& a : va) {
          ARTS_USER_ERROR_IF(
              a.get_name() not_eq va.front().get_name(),
              "An ArrayOfAgenda must only consist of agendas with the same name\n"
              "You have input a list of agendas that contains disimilar names.\n"
              "\nThe first item is named: \"",
              va.front().get_name(),
              '"',
              '\n',
              "A later item in the list is names: \"",
              a.get_name(),
              '"',
              '\n')
        }
        return va;
      }), "Create from :class:`list`")
      .def("__repr__", [](ArrayOfAgenda&) { return "ArrayOfAgenda"; })
      .def("__str__",
           [](ArrayOfAgenda& aa) {
             return var_string(aa);
           })
      .PythonInterfaceFileIO(ArrayOfAgenda)
      .def("__len__", [](const ArrayOfAgenda& x) { return x.nelem(); })
      .def(
          "__getitem__",
          [](ArrayOfAgenda& x, Index i) -> Agenda& {
            if (x.nelem() <= i or i < 0)
              throw std::out_of_range(var_string("Bad index access: ",
                                                 i,
                                                 " in object of size [0, ",
                                                 x.size(),
                                                 ")"));
            return x[i];
          },
          py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](ArrayOfAgenda& x, Index i, Agenda y) {
            if (x.nelem() <= i or i < 0) {
              throw std::out_of_range(var_string("Bad index access: ",
                                                 i,
                                                 " in object of size [0, ",
                                                 x.size(),
                                                 ")"));
            }
            if (y.get_name() not_eq x.front().get_name()) y.set_name(x.front().get_name());
            x[i] = std::move(y);
          },
          py::return_value_policy::reference_internal)
      .def("append",
           [](ArrayOfAgenda& x, Agenda y) {
             if (x.nelem() and y.get_name() not_eq x.front().get_name())
               y.set_name(x.front().get_name());
             x.emplace_back(std::move(y));
           }, "Appends a :class:`~pyarts.arts.Agenda` at the end of the array")
      .def(
          "pop",
          [](ArrayOfAgenda& aa) {
            Agenda a = std::move(aa.back());
            aa.pop_back();
            return a;
          },
          "Pops and returns the last item")
      .def(
          "finalize",
          [](ArrayOfAgenda& aa) {
            for (auto& a : aa)
              a.finalize();
          },
          py::doc("Checks if the agenda works"))
      .def_property(
          "name",
          [](ArrayOfAgenda& a) -> String {
            if (a.nelem() == 0) return "";
            return a.front().get_name();
          },
          [](ArrayOfAgenda& aa, const String& name) {
            for (auto& a : aa) a.set_name(name);
          }, py::doc(":class:`~pyarts.arts.String` Name of the array of agenda"))
      .def(py::pickle(
          [](const ArrayOfAgenda& v) {
            auto n = v.size();
            std::vector<Agenda> out(n);
            std::copy(v.begin(), v.end(), out.begin());
            return py::make_tuple(std::move(out));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_shared<ArrayOfAgenda>(t[0].cast<std::vector<Agenda>>());
          }))
      .PythonInterfaceWorkspaceDocumentationExtra(ArrayOfAgenda, "\n\nThese arrays are partial to contain inter-item logic, please be cautious using them");
  py::implicitly_convertible<std::vector<Agenda>, ArrayOfAgenda>();
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize agendas\n", e.what()));
}
}  // namespace Python
