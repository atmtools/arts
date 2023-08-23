#include <parameters.h>
#include <python_interface.h>

#include <algorithm>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "auto_wsv.h"
#include "workspace_agenda_class.h"
#include "workspace_class.h"

extern Parameters parameters;
namespace Python {
py::tuple pickle_method(const Method& m) {
  return py::make_tuple(
      m.get_name(), m.get_ins(), m.get_outs(), m.get_setval(), m.overwrite());
}

Method unpickle_method(const py::tuple& t) {
  ARTS_USER_ERROR_IF(t.size() != 5, "Invalid state!")

  auto name = t[0].cast<std::string>();
  auto ins = t[1].cast<std::vector<std::string>>();
  auto outs = t[2].cast<std::vector<std::string>>();
  auto setval = t[3].cast<std::optional<Wsv>>();
  auto ow = t[4].cast<bool>();

  return Method{name, ins, outs, setval, ow};
}

py::tuple pickle_agenda(const Agenda& ag) {
  return py::make_tuple(ag.get_name(),
                        ag.get_methods(),
                        ag.get_share(),
                        ag.get_copy(),
                        ag.is_checked());
}

Agenda unpickle_agenda(const py::tuple& t) {
  ARTS_USER_ERROR_IF(t.size() != 5, "Invalid state!")

  auto n = t[0].cast<std::string>();
  auto m = t[1].cast<std::vector<Method>>();
  auto s = t[2].cast<std::vector<std::string>>();
  auto c = t[3].cast<std::vector<std::string>>();
  auto ch = t[4].cast<bool>();

  return Agenda{n, m, s, c, ch};
}

std::filesystem::path correct_include_path(
    const std::filesystem::path& path_copy);

void py_auto_wsv(artsclass<Workspace>& ws);
void py_auto_wsm(artsclass<Workspace>& ws);

void py_workspace(artsclass<Workspace>& ws) try {
  ws.def(py::init([](bool with_defaults) {
           if (with_defaults)
             return std::make_shared<Workspace>(WorkspaceInitialization::FromGlobalDefaults);
           return std::make_shared<Workspace>(WorkspaceInitialization::Empty);
         }),
         py::arg("with_defaults") = true)
      .def(py::init([](Workspace& w) { return w; }))
      .def("__copy__", [](Workspace& w) { return w; })
      .def("__deepcopy__", [](Workspace& w, py::dict&) { return w; })
      .def(
          "get",
          [](Workspace& w, const std::string& n) { return from(w.share(n)); },
          py::return_value_policy::reference_internal,
          py::keep_alive<0, 1>())
      .def("set",
           [](Workspace& w, const std::string& n, const PyWsvValue& x) {
             w.set(n, std::make_shared<Wsv>(from(x)));

             auto ptr = w.share(n);
             if (ptr->holds<Agenda>()) {
               auto& ag = ptr->get_unsafe<Agenda>();
               if (not ag.is_checked()) {
                 ag.set_name(n);
                 ag.finalize();
               }
             } else if (ptr->holds<ArrayOfAgenda>()) {
               auto& ags = ptr->get_unsafe<ArrayOfAgenda>();
               for (auto& ag : ags) {
                if (not ag.is_checked()) {
                  ag.set_name(n);
                  ag.finalize();
                }
               }
             }
           })
      .def("has",
           [](Workspace& w, const std::string& n) { return w.contains(n); });

  ws.def("__str__", [](const Workspace& w) {
    return var_string(w);
  });

  py_auto_wsv(ws);
  py_auto_wsm(ws);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize workspace\n", e.what()));
}
}  // namespace Python
