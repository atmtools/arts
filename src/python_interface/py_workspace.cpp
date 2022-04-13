#include <global_data.h>
#include <parameters.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl/filesystem.h>

#include "py_auto_interface.h"
#include "py_macros.h"
#include "python_interface.h"

extern Parameters parameters;

void parse_path_from_environment(String envvar, ArrayOfString& paths);
namespace Python {
std::filesystem::path correct_include_path(
    const std::filesystem::path& path_copy);
Agenda* parse_agenda(const char* filename, const Verbosity& verbosity);
void py_auto_workspace(py::class_<Workspace>&, py::class_<WorkspaceVariable>&);

void py_workspace(py::module_& m,
                  py::class_<Workspace>& ws,
                  py::class_<WorkspaceVariable>& wsv) {
  static bool init = true;
  if (init) {
    init = false;

    define_wsv_group_names();
    Workspace::define_wsv_data();
    Workspace::define_wsv_map();
    define_md_data_raw();
    expand_md_data_raw_to_md_data();
    define_md_map();
    define_md_raw_map();
    define_agenda_data();
    define_agenda_map();
    ARTS_ASSERT(check_agenda_data());
    global_data::workspace_memory_handler.initialize();

    // Set parameters that are know on first execution
#ifdef ARTS_DEFAULT_INCLUDE_DIR
    String arts_default_include_path(ARTS_DEFAULT_INCLUDE_DIR);
    if (arts_default_include_path != "" && !parameters.includepath.nelem()) {
      // Skip delimiters at beginning.
      String::size_type lastPos =
          arts_default_include_path.find_first_not_of(":", 0);
      // Find first "non-delimiter".
      String::size_type pos =
          arts_default_include_path.find_first_of(":", lastPos);

      while (String::npos != pos || String::npos != lastPos) {
        parameters.includepath.push_back(
            arts_default_include_path.substr(lastPos, pos - lastPos));
        lastPos = arts_default_include_path.find_first_not_of(":", pos);
        pos = arts_default_include_path.find_first_of(":", lastPos);
      }
    }
#endif

    parse_path_from_environment("ARTS_INCLUDE_PATH", parameters.includepath);
    parse_path_from_environment("ARTS_DATA_PATH", parameters.datapath);

    parameters.includepath.insert(parameters.includepath.begin(), ".");
    parameters.datapath.insert(parameters.datapath.begin(), ".");
  }

  ws.def(py::init([](Index verbosity, Index agenda_verbosity) {
           auto* w = new Workspace{};
           w->initialize();
           w->push(w->WsvMap.at("verbosity"),
                   new Verbosity{agenda_verbosity, verbosity, 0});
           return w;
         }),
         py::arg("verbosity") = 0,
         py::arg("agenda_verbosity") = 0)
      .def("execute_controlfile",
           [](Workspace& w, const std::filesystem::path& path) {
             std::unique_ptr<Agenda> a{parse_agenda(
                 correct_include_path(path).c_str(),
                 *reinterpret_cast<Verbosity*>(w[w.WsvMap.at("verbosity")]))};
             w.initialize();
             a->execute(w);
           })
      .def_property_readonly("size", &Workspace::nelem)
      .def(
          "create_variable",
          [](Workspace& w,
             const char* group,
             const char* name_,
             std::optional<const char*> desc) {
            using global_data::workspace_memory_handler;

            const Index group_index = get_wsv_group_id(group);
            ARTS_USER_ERROR_IF(
                group_index < 0, "Cannot recognize group: ", group)

            String name;
            if (name_) {
              name = name_;
            } else {
              Index i = 0;
              while (w.WsvMap.find(
                         name = var_string("::auto::", group, i++)) not_eq
                     w.WsvMap.end()) {
              }
            }

            auto varptr = w.WsvMap.find(name);
            ARTS_USER_ERROR_IF(
                varptr not_eq w.WsvMap.end() and
                    w.wsv_data.at(varptr->second).Group() not_eq group_index,
                "Already exist of different group: ",
                name)

            const Index pos = w.add_wsv_inplace(
                WsvRecord(name.c_str(),
                          desc.has_value() ? *desc : "No description",
                          group_index));
            w.push(pos, workspace_memory_handler.allocate(group_index));
            return WorkspaceVariable{w, pos};
          },
          py::keep_alive<0, 1>(),
          py::arg("group").none(false),
          py::arg("name"),
          py::arg("desc") = std::nullopt);

  ws.def("_hasattr_check_", [](Workspace& w, const char* name) -> bool {
    return w.WsvMap.find(name) != w.WsvMap.end();
  });

  ws.def(
      "_getattr_unchecked_",
      [](Workspace& w, const char* name) -> WorkspaceVariable {
        return WorkspaceVariable{w, w.WsvMap.at(name)};
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
        if (auto x = w.wsv_data.at(i).Name(); x[0] != ':')
          vars.push_back(var_string(
              x,
              ": ",
              global_data::wsv_group_names.at(w.wsv_data.at(i).Group())));
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
      .def_property_readonly("name", &WorkspaceVariable::name)
      .def_property_readonly("desc",
                             [](WorkspaceVariable& wv) {
                               return wv.ws.wsv_data.at(wv.pos).Description();
                             })
      .def_property_readonly("init", &WorkspaceVariable::is_initialized)
      .def("initialize_if_not", &WorkspaceVariable::initialize_if_not)
      .def("__str__",
           [](const WorkspaceVariable& w) {
             return var_string("Workspace ",
                               global_data::wsv_group_names.at(
                                   w.ws.wsv_data.at(w.pos).Group()));
           })
      .def("__repr__", [](const WorkspaceVariable& w) {
        return var_string(
            "Workspace ",
            global_data::wsv_group_names.at(w.ws.wsv_data.at(w.pos).Group()));
      });

  py::class_<WsvRecord>(m, "WsvRecord")
      .PythonInterfaceCopyValue(WsvRecord)
      .def_property_readonly("name", &WsvRecord::Name)
      .def_property_readonly("description", &WsvRecord::Description)
      .def_property_readonly("group", &WsvRecord::Group)
      .def_property_readonly("implicit", &WsvRecord::Implicit)
      .def("__str__", [](const WsvRecord& mr) { return var_string(mr); })
      .def("__repr__", [](const WsvRecord& mr) { return var_string(mr); });

  py::class_<Array<WsvRecord>>(m, "ArrayOfWsvRecord")
      .PythonInterfaceArrayDefault(WsvRecord)
      .PythonInterfaceBasicRepresentation(Array<WsvRecord>);

  // Should be last as it contains several implicit conversions
  py_auto_workspace(ws, wsv);
}
}  // namespace Python