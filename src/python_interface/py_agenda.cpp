#include <algorithm>
#include <parameters.h>
#include <pybind11/cast.h>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>

#include "auto_wsg.h"
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
  artsclass<CallbackOperator>(m, "CallbackOperator")
      .def(py::init([](const std::function<void(Workspace&)>& f,
                       const std::vector<std::string>& i,
                       const std::vector<std::string>& o) {
             return std::make_shared<CallbackOperator>(f, i, o);
           }),
           py::arg("f") = std::function<void(Workspace&)>([](Workspace&) { throw std::runtime_error("No-op"); }),
           py::arg("inputs")=std::vector<std::string>{},
           py::arg("outputs")=std::vector<std::string>{},
           py::doc("Initialize as structured call"))
      .PythonInterfaceCopyValue(CallbackOperator)
      .PythonInterfaceBasicRepresentation(CallbackOperator)
      .def("__call__", [](CallbackOperator& f, Workspace& ws) { f(ws); })
      .PythonInterfaceWorkspaceDocumentation(CallbackOperator);

  artsclass<Method>(m, "Method")
      .def(py::init([](const std::string& n,
                       const std::vector<std::string>& a,
                       const std::unordered_map<std::string, std::string>& kw)
                        -> Method {
             return {n, a, kw};
           }),
           py::arg("name"),
           py::arg("args") = std::vector<std::string>{},
           py::arg("kwargs") = std::unordered_map<std::string, std::string>{},
           py::doc("A named method with args and kwargs"))
      .def(py::init([](const std::string& n, const PyWsvValue& v) -> Method {
             return std::visit(
                 [&](auto&& wsv) -> Method {
                   return {n, Wsv{*wsv}, n.front() == '@'};
                 },
                 from(v).value);
           }),
           py::arg("name"),
           py::arg("wsv"),
           py::doc("A method that sets a workspace variable"))
      .def_property_readonly(
          "val",
          [](const Method& method) -> std::variant<py::none, PyWsvValue> {
            const auto& x = method.get_setval();
            if (x) return from(x.value());
            return py::none();
          },
          py::doc("The value (if any) of a set method"))
      .def_property_readonly(
          "name",
          [](const Method& method) { return method.get_name(); },
          py::doc("The name of the method"))
      .def("__str__", [](const Method& method){return var_string(method);})
      .doc() = "The method class of ARTS";

  artsclass<Agenda>(m, "Agenda")
      .def(py::init([]() {
    return std::make_shared<Agenda>(); }), "Create empty")
      .def(py::init<std::string>(), py::arg("name"), "Create with name")
      .PythonInterfaceWorkspaceVariableConversion(Agenda)
      .PythonInterfaceFileIO(Agenda)
      .def("__copy__", [](const Agenda& x) -> std::shared_ptr<Agenda>{return std::make_shared<Agenda>(x);})
      .def("__deepcopy__", [](const Agenda& x, py::dict&) -> std::shared_ptr<Agenda>{return std::make_shared<Agenda>(x);})
      .def(
          "add",
          &Agenda::add,
          py::arg("method").none(false),
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
          [](Agenda& a, Workspace& ws) {
            a.execute(ws);
          }, py::arg("ws"),
          "Executes the agenda on the provided workspace")
      .def(
          "finalize",
          [](Agenda& a, bool fix) {
            a.finalize(fix);
          }, py::arg("fix") = false, "Finalize the agenda, making it possible to use it in the workspace")
      .def_property_readonly(
          "name",
          [](const Agenda& agenda) { return agenda.get_name(); },
          py::doc("The name of the agenda"))
      .def_property_readonly(
          "methods",
          [](const Agenda& agenda) { return agenda.get_methods(); },
          py::doc("The methods of the agenda"))
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
        return std::make_shared<ArrayOfAgenda>(va);
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
