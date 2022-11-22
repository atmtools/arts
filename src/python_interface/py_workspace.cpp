#include <global_data.h>
#include <memory>
#include <mutex>
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
#include "tokval.h"

extern Parameters parameters;
namespace Python {
py::tuple pickle_agenda(const Agenda &ag);
Agenda unpickle_agenda(Workspace &ws, const py::tuple &t);

py::tuple pickle_mrecord(const MRecord &mr) {
  return py::make_tuple(
      mr.Id(), pickle_agenda(mr.Tasks()), mr.SetValue(),
      mr.Tasks().workspace()->wsvs(mr.Out()),
      mr.Tasks().workspace()->wsvs(mr.In()),
      global_data::md_data[mr.Id()].Name(),
      global_data::md_data[mr.Id()].Description(),
      global_data::md_data[mr.Id()].Authors(),
      global_data::md_data[mr.Id()].Out(), global_data::md_data[mr.Id()].GOut(),
      global_data::md_data[mr.Id()].GOutType(),
      global_data::md_data[mr.Id()].GOutSpecType(),
      global_data::md_data[mr.Id()].GOutDescription(),
      global_data::md_data[mr.Id()].In(), global_data::md_data[mr.Id()].GIn(),
      global_data::md_data[mr.Id()].GInType(),
      global_data::md_data[mr.Id()].GInSpecType(),
      global_data::md_data[mr.Id()].GInDefault(),
      global_data::md_data[mr.Id()].GInDescription(),
      global_data::md_data[mr.Id()].InOnly(),
      global_data::md_data[mr.Id()].InOut(),
      global_data::md_data[mr.Id()].OutOnly(),
      global_data::md_data[mr.Id()].SetMethod(),
      global_data::md_data[mr.Id()].AgendaMethod(),
      global_data::md_data[mr.Id()].Supergeneric(),
      global_data::md_data[mr.Id()].UsesTemplates(),
      global_data::md_data[mr.Id()].PassWorkspace(),
      global_data::md_data[mr.Id()].PassWsvNames(),
      global_data::md_data[mr.Id()].ActualGroups());
}

MRecord unpickle_mrecord(Workspace &ws, const py::tuple &t) {
  ARTS_USER_ERROR_IF(t.size() != 29, "Invalid state!")

  auto id = t[0].cast<Index>();
  auto ag_tup = t[1].cast<py::tuple>();
  auto tv = t[2].cast<TokVal>();
  auto out_wsv = t[3].cast<Workspace::wsv_data_type>();
  auto in_wsv = t[4].cast<Workspace::wsv_data_type>();

  ARTS_USER_ERROR_IF(global_data::md_data[id].Name() not_eq t[5].cast<String>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].Description() not_eq
                         t[6].cast<String>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].Authors() not_eq
                         t[7].cast<ArrayOfString>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].Out() not_eq
                         t[8].cast<ArrayOfIndex>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].GOut() not_eq
                         t[9].cast<ArrayOfString>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].GOutType() not_eq
                         t[10].cast<ArrayOfIndex>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].GOutSpecType() not_eq
                         t[11].cast<ArrayOfArrayOfIndex>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].GOutDescription() not_eq
                         t[12].cast<Array<String>>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].In() not_eq
                         t[13].cast<ArrayOfIndex>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].GIn() not_eq
                         t[14].cast<ArrayOfString>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].GInType() not_eq
                         t[15].cast<ArrayOfIndex>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].GInSpecType() not_eq
                         t[16].cast<ArrayOfArrayOfIndex>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].GInDefault() not_eq
                         t[17].cast<Array<String>>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].GInDescription() not_eq
                         t[18].cast<Array<String>>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].InOnly() not_eq
                         t[19].cast<ArrayOfIndex>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].InOut() not_eq
                         t[20].cast<ArrayOfIndex>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].OutOnly() not_eq
                         t[21].cast<ArrayOfIndex>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].SetMethod() not_eq
                         t[22].cast<bool>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].AgendaMethod() not_eq
                         t[23].cast<bool>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].Supergeneric() not_eq
                         t[24].cast<bool>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].UsesTemplates() not_eq
                         t[25].cast<bool>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].PassWorkspace() not_eq
                         t[26].cast<bool>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].PassWsvNames() not_eq
                         t[27].cast<bool>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")
  ARTS_USER_ERROR_IF(global_data::md_data[id].ActualGroups() not_eq
                         t[28].cast<String>(),
                     "Method state not same as on construction; are you using "
                     "a different Arts version?")

  return MRecord{id, ws.wsvs(out_wsv), ws.wsvs(in_wsv), tv,
                 unpickle_agenda(ws, ag_tup)};
}

py::tuple pickle_agenda(const Agenda &ag) {
  std::vector<py::tuple> methods(ag.Methods().size());
  std::transform(ag.Methods().begin(), ag.Methods().end(), methods.begin(),
                 [](auto &in) { return pickle_mrecord(in); });
  return py::make_tuple(ag.name(), methods, ag.is_main_agenda());
}

Agenda unpickle_agenda(Workspace &ws, const py::tuple &t) {
  ARTS_USER_ERROR_IF(t.size() != 3, "Invalid state!")

  auto val = Agenda{ws};

  val.set_name(t[0].cast<String>());
  auto mr_tup = t[1].cast<std::vector<py::tuple>>();
  Array<MRecord> methods(mr_tup.size());
  std::transform(mr_tup.begin(), mr_tup.end(), methods.begin(),
                 [&ws](auto &met) { return unpickle_mrecord(ws, met); });
  val.set_methods(methods);
  val.set_outputs_to_push_and_dup(Verbosity{});

  if (t[2].cast<bool>())
    val.set_main_agenda();
  else
    val.check(ws, Verbosity{});

  return val;
}

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
           auto w = Workspace::create();
           w->push_move(w->WsvMap_ptr->at("verbosity"),
                   std::make_shared<Verbosity>(agenda_verbosity, verbosity, 0));
           return w;
         }),
         py::arg("verbosity") = 0,
         py::arg("agenda_verbosity") = 0)
      .def(py::init([](Workspace& w) {return w.shallowcopy();}))
      .def(
          "__copy__",
          [](Workspace& w) { return w.shallowcopy(); },
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
        std::map<String, py::object> value;
        std::map<String, py::tuple> agenda;
        std::map<String, std::vector<py::tuple>> array_agenda;

        for (Index i = 0; i < w.nelem(); i++) {
          auto &wsv_data = w.wsv_data_ptr->operator[](i);

          ARTS_USER_ERROR_IF(
              wsv_data.Group() == WorkspaceGroupIndexValue<CallbackFunction>,
              "Cannot pickle a workspace that contains a CallbackFunction")
          if (w.is_initialized(i)) {
            const String &name = wsv_data.Name();

            if (wsv_data.Group() == WorkspaceGroupIndexValue<Agenda>) {
              agenda[name] = pickle_agenda(WorkspaceVariable{w, i});
            } else if (wsv_data.Group() ==
                       WorkspaceGroupIndexValue<ArrayOfAgenda>) {
              const ArrayOfAgenda &ags = WorkspaceVariable{w, i};
              std::vector<py::tuple> &ags_tup = array_agenda[name];
              ags_tup.resize(ags.size());
              std::transform(ags.begin(), ags.end(), ags_tup.begin(),
                             [](auto &in) { return pickle_agenda(in); });
            } else {
              value[name] =
                  py::type::of<WorkspaceVariable>()(WorkspaceVariable{w, i})
                      .attr("value");
            }
          }
        }

        return py::make_tuple(value,
                              agenda,
                              array_agenda,
                              *w.wsv_data_ptr,
                              global_data::AgendaMap,
                              global_data::MdMap,
                              global_data::WsvMap,
                              global_data::WsvGroupMap,
                              get_group_types());
      },
      [](const py::tuple& t) {
        auto out = Workspace::create();

        ARTS_USER_ERROR_IF(t.size() != 9, "Invalid state!")

        ARTS_USER_ERROR_IF((t[4].cast<std::map<String, Index>>() not_eq
                            global_data::AgendaMap),
                           "Mismatch between Agendas"
                           " at time of pickling and unpickling the workspace")
        ARTS_USER_ERROR_IF(
            (t[5].cast<std::map<String, Index>>() not_eq global_data::MdMap),
            "Mismatch between Methods"
            " at time of pickling and unpickling the workspace")
        ARTS_USER_ERROR_IF(
            (t[6].cast<std::map<String, Index>>() not_eq global_data::WsvMap),
            "Mismatch between Workspace Variable Names"
            " at time of pickling and unpickling the workspace")
        ARTS_USER_ERROR_IF((t[7].cast<std::map<String, Index>>() not_eq
                            global_data::WsvGroupMap),
                           "Mismatch between Groups"
                           " at time of pickling and unpickling the workspace")
        ARTS_USER_ERROR_IF((t[8].cast<ArrayOfIndex>() not_eq get_group_types()),
                           "Mismatch between Workspace Variable Groups"
                           " at time of pickling and unpickling the workspace")

        auto value = t[0].cast<std::map<String, py::object>>();
        auto agenda = t[1].cast<std::map<String, py::tuple>>();
        auto array_agenda = t[2].cast<std::map<String, std::vector<py::tuple>>>();
        out->wsvs(t[3].cast<Workspace::wsv_data_type>());

        for (auto &x : value) {
          auto wsv_ptr = out->WsvMap_ptr->find(x.first);
          ARTS_USER_ERROR_IF(wsv_ptr == out->WsvMap_ptr->end(),
                             "Cannot find WSV ", x.first)
          py::type::of<WorkspaceVariable>()(
              WorkspaceVariable{*out, wsv_ptr->second})
              .attr("value") = x.second;
        }

        for (auto &x : agenda) {
          auto wsv_ptr = out->WsvMap_ptr->find(x.first);
          ARTS_USER_ERROR_IF(wsv_ptr == out->WsvMap_ptr->end(),
                             "Cannot find WSV ", x.first)
          py::type::of<WorkspaceVariable>()(
              WorkspaceVariable{*out, wsv_ptr->second})
              .attr("value") = unpickle_agenda(*out, x.second);
        }

        for (auto &x : array_agenda) {
          auto wsv_ptr = out->WsvMap_ptr->find(x.first);
          ARTS_USER_ERROR_IF(wsv_ptr == out->WsvMap_ptr->end(),
                             "Cannot find WSV ", x.first)
          ArrayOfAgenda agendas(x.second.size());
          std::transform(
              x.second.begin(), x.second.end(), agendas.begin(),
              [&out](auto &in) { return unpickle_agenda(*out, in); });
          py::type::of<WorkspaceVariable>()(
              WorkspaceVariable{*out, wsv_ptr->second})
              .attr("value") = agendas;
        }

        return out;
      }));

  // Should be last as it contains several implicit conversions
  py_auto_workspace(ws, wsv);
}
}  // namespace Python
