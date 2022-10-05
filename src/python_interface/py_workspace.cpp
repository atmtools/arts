#include <global_data.h>
#include <parameters.h>
#include <wsv_aux_operator.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl/filesystem.h>
#include <algorithm>
#include <vector>

#include "debug.h"
#include "py_auto_interface.h"
#include "py_macros.h"
#include "python_interface.h"

extern Parameters parameters;
namespace Python {
std::filesystem::path correct_include_path(
    const std::filesystem::path& path_copy);

Agenda* parse_agenda(Workspace&,
                     const char* filename,
                     const Verbosity& verbosity);

void py_auto_workspace(py::class_<Workspace, std::shared_ptr<Workspace>>&,
                       py::class_<WorkspaceVariable>&);

auto get_group_types() {
  ArrayOfIndex tmp(global_data::wsv_data.nelem());
  std::transform(global_data::wsv_data.begin(),
                 global_data::wsv_data.end(),
                 tmp.begin(),
                 [](auto& wsv) { return wsv.Group(); });
  return tmp;
};

void py_workspace(py::module_& m,
                  py::class_<Workspace, std::shared_ptr<Workspace>>& ws,
                  py::class_<WorkspaceVariable>& wsv) {
  ws.def(py::init([](Index verbosity, Index agenda_verbosity) {
           auto* w = new Workspace{};
           w->push_move(w->WsvMap_ptr->at("verbosity"),
                   std::make_shared<Verbosity>(agenda_verbosity, verbosity, 0));
           return w;
         }),
         py::arg("verbosity") = 0,
         py::arg("agenda_verbosity") = 0)
      .def(py::init([](Workspace& w) {return new Workspace{w};}))
      .def(
          "__copy__",
          [](Workspace& w) -> Workspace { return w; },
          py::is_operator())
      .def(
          "__deepcopy__",
          [](Workspace& w, py::dict&) {
            return w.deepcopy();
          },
          py::is_operator())
      .def("execute_controlfile",
           [](Workspace& w, const std::filesystem::path& path) {
             std::unique_ptr<Agenda> a{parse_agenda(w, 
                 correct_include_path(path).c_str(),
                 *static_cast<Verbosity*>(w.get<Verbosity>("verbosity")))};
             a->execute(w);
           })
      .def_property_readonly("size", &Workspace::nelem)
      .def(
          "create_variable",
          [](Workspace& w,
             const char* group,
             const char* name_,
             std::optional<const char*> desc) {
            const Index group_index = get_wsv_group_id(group);
            ARTS_USER_ERROR_IF(
                group_index < 0, "Cannot recognize group: ", group)

            String name;
            if (name_) {
              name = name_;
            } else {
              Index i = 0;
              while (w.WsvMap_ptr->find(
                         name = var_string("::auto::", group, i++)) not_eq
                     w.WsvMap_ptr->end()) {
              }
            }

            auto varptr = w.WsvMap_ptr->find(name);
            ARTS_USER_ERROR_IF(
                varptr not_eq w.WsvMap_ptr->end() and
                    w.wsv_data_ptr->at(varptr->second).Group() not_eq group_index,
                "Already exist of different group: ",
                name)

            const Index pos = w.add_wsv(
                WsvRecord(name.c_str(),
                          desc.has_value() ? *desc : "No description",
                          group_index));
            w.emplace(pos);
            return WorkspaceVariable{w, pos};
          },
          py::keep_alive<0, 1>(),
          py::arg("group").none(false),
          py::arg("name"),
          py::arg("desc") = std::nullopt);

  ws.def("_hasattr_check_", [](Workspace& w, const char* name) -> bool {
    return w.WsvMap_ptr->find(name) != w.WsvMap_ptr->end();
  });

  ws.def(
      "_getattr_unchecked_",
      [](Workspace& w, const char* name) -> WorkspaceVariable {
        return WorkspaceVariable{w, w.WsvMap_ptr->at(name)};
      },
      py::keep_alive<0, 1>());

  ws.def(
      "swap",
      [](Workspace& w_1, Workspace& w_2) { w_1.swap(w_2); },
      py::arg("ws").noconvert().none(false));

  ws.def("__repr__", [](const Workspace&) { return "Workspace"; });

  ws.def("__str__", [](const Workspace& w) {
    ArrayOfString vars;
    for (Index i = 0; i < w.nelem(); i++)
      if (w.is_initialized(i))
        if (auto x = w.wsv_data_ptr->at(i).Name(); x[0] != ':')
          vars.push_back(var_string(
              x,
              ": ",
              global_data::wsv_groups.at(w.wsv_data_ptr->at(i).Group())));
    std::sort(vars.begin(), vars.end());
    return var_string("Workspace [ ", stringify(vars, ", "), ']');
  });

  ws.PythonInterfaceCopyValue(Workspace);

  wsv.def(py::init([](WorkspaceVariable& w) { return w; }))
      .def(
          "__copy__",
          [](WorkspaceVariable&) { ARTS_USER_ERROR("Cannot copy") },
          "Cannot copy a workspace variable, try instead to copy its value")
      .def(
          "__deepcopy__",
          [](WorkspaceVariable&, py::dict&) { ARTS_USER_ERROR("Cannot copy") },
          "Cannot copy a workspace variable, try instead to copy its value")
      .def_property_readonly("name", &WorkspaceVariable::name, "The name of the workspace variable")
      .def_property_readonly("desc",
                             [](WorkspaceVariable& wv) {
                               return wv.ws.wsv_data_ptr->at(wv.pos).Description();
                             }, "A description of the workspace variable")
      .def_property_readonly("init", &WorkspaceVariable::is_initialized, "Flag if the variable is initialized")
      .def("initialize_if_not", &WorkspaceVariable::initialize_if_not, "Default initialize the variable if it is not already initialized")
      .def("__str__",
           [](const WorkspaceVariable& w) {
             return var_string("Workspace ",
                               global_data::wsv_groups.at(w.group()));
           })
      .def("__repr__",
           [](const WorkspaceVariable& w) {
             return var_string("Workspace ",
                               global_data::wsv_groups.at(w.group()));
           })
      .def_property_readonly("group",
                             [](const WorkspaceVariable& w) {
                               return global_data::wsv_groups.at(w.group());
                             }, "The group of the variable")
      .def("delete_level", [](WorkspaceVariable& v) {
        if (v.stack_depth()) v.pop_workspace_level();
      }, "Delete a level from the workspace variable stack");

  ws.def("number_of_initialized_variables", [](Workspace& w){
    Index count = 0;
    for (Index i=0; i<w.nelem(); i++) count += w.is_initialized(i);
    return count;
  });

  py::class_<WsvRecord>(m, "WsvRecord")
      .PythonInterfaceCopyValue(WsvRecord)
      .def_property_readonly("name", &WsvRecord::Name)
      .def_property_readonly("description", &WsvRecord::Description)
      .def_property_readonly("group", &WsvRecord::Group)
      .def_property_readonly("defval", &WsvRecord::default_value)
      .def_property_readonly("groupname",
                             [](const WsvRecord& x) {
                               return global_data::wsv_groups[x.Group()].name;
                             })
      .def("__str__", [](const WsvRecord& mr) { return var_string(mr); })
      .def("__repr__", [](const WsvRecord& mr) { return var_string(mr); })
      .def(py::pickle(
          [](const WsvRecord& self) {
            return py::make_tuple(self.Name(),
                                  self.Description(),
                                  global_data::wsv_groups[self.Group()].name,
                                  self.default_value());
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 4, "Invalid state!")
            return new WsvRecord{t[0].cast<String>().c_str(),
                                 t[1].cast<String>().c_str(),
                                 t[2].cast<String>(),
                                 t[3].cast<TokVal>()};
          }));

  py::class_<Array<WsvRecord>>(m, "ArrayOfWsvRecord")
      .PythonInterfaceArrayDefault(WsvRecord)
      .PythonInterfaceBasicRepresentation(Array<WsvRecord>);

  ws.def(py::pickle(
      [](Workspace& w) {
        std::vector<py::object> value;
        value.reserve(w.nelem());

        for (Index i = 0; i < w.nelem(); i++) {
          auto& wsv_data = w.wsv_data_ptr->operator[](i);

          ARTS_USER_ERROR_IF(
              wsv_data.Group() == WorkspaceGroupIndexValue<CallbackFunction>,
              "Cannot pickle a workspace that contains a CallbackFunction")

          if (not w.is_initialized(i)) {
            value.push_back(py::none{});
          } else {
            value.push_back(
                py::type::of<WorkspaceVariable>()(WorkspaceVariable{w, i})
                    .attr("value"));
          }
        }

        return py::make_tuple(value,
                              *w.wsv_data_ptr,
                              global_data::AgendaMap,
                              global_data::MdMap,
                              global_data::WsvMap,
                              global_data::WsvGroupMap,
                              get_group_types());
      },
      [](const py::tuple& t) {
        ARTS_USER_ERROR_IF(t.size() != 7, "Invalid state!")

        ARTS_USER_ERROR_IF((t[2].cast<std::map<String, Index>>() not_eq
                            global_data::AgendaMap),
                           "Mismatch between Agendas"
                           " at time of pickling and unpickling the workspace")
        ARTS_USER_ERROR_IF(
            (t[3].cast<std::map<String, Index>>() not_eq global_data::MdMap),
            "Mismatch between Methods"
            " at time of pickling and unpickling the workspace")
        ARTS_USER_ERROR_IF(
            (t[4].cast<std::map<String, Index>>() not_eq global_data::WsvMap),
            "Mismatch between Workspace Variable Names"
            " at time of pickling and unpickling the workspace")
        ARTS_USER_ERROR_IF((t[5].cast<std::map<String, Index>>() not_eq
                            global_data::WsvGroupMap),
                           "Mismatch between Groups"
                           " at time of pickling and unpickling the workspace")
        ARTS_USER_ERROR_IF((t[6].cast<ArrayOfIndex>() not_eq get_group_types()),
                           "Mismatch between Workspace Variable Groups"
                           " at time of pickling and unpickling the workspace")

        auto out = std::shared_ptr<Workspace>(new Workspace{});

        const auto aoi = out->wsvs(t[1].cast<Workspace::wsv_data_type>());
        auto value = t[0].cast<std::vector<py::object>>();

        for (auto i : aoi) {
          if (py::isinstance<Agenda>(value[i])) continue;
          if (py::isinstance<ArrayOfAgenda>(value[i])) continue;
          if (py::isinstance<py::none>(value[i])) continue;
          py::type::of<WorkspaceVariable>()(WorkspaceVariable{*out, i})
              .attr("value") = value[i];
        }

        for (auto i : aoi) {
          if (not py::isinstance<Agenda>(value[i]) or
              py::isinstance<py::none>(value[i]))
            continue;
          py::type::of<WorkspaceVariable>()(WorkspaceVariable{*out, i})
              .attr("value") = value[i].cast<Agenda>().deepcopy_if(*out);
        }

        for (auto i : aoi) {
          if (not py::isinstance<ArrayOfAgenda>(value[i]) or
              py::isinstance<py::none>(value[i]))
            continue;
          py::type::of<WorkspaceVariable>()(WorkspaceVariable{*out, i})
              .attr("value") =
              deepcopy_if(*out, value[i].cast<ArrayOfAgenda>());
        }

        return out;
      }));

  // Should be last as it contains several implicit conversions
  py_auto_workspace(ws, wsv);
}
}  // namespace Python