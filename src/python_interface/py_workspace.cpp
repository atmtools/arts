
#include <py_auto_interface.h>

#include <global_data.h>
#include <parameters.h>
#include <workspace_ng.h>

#include "py_macros.h"

#include <cstdlib>

extern Parameters parameters;
void parse_path_from_environment(String envvar, ArrayOfString &paths);

namespace Python {

void py_auto_workspace(py::class_<Workspace>&);
void py_workspace(py::module_& m) {
  static bool init=true;
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

  auto ws = py::class_<Workspace>(m, "Workspace").def(py::init([]() {
    Workspace w{};
    w.initialize();
    w.push(Workspace::WsvMap.at("verbosity"), new Verbosity{});
    return w;
  }));

  ws.def_property_readonly("size", &Workspace::nelem);

  ws.def(
      "create_variable",
      [](Workspace& w, char* group, char* name, char* desc) {
        using global_data::workspace_memory_handler;

        ARTS_USER_ERROR_IF(std::find_if(Workspace::WsvMap.begin(),
                                        Workspace::WsvMap.end(),
                                        [name](auto& b) {
                                          return b.first == name;
                                        }) not_eq Workspace::WsvMap.end(),
                           "A variable of this name already exist: ",
                           name)

        const Index group_index = get_wsv_group_id(group);
        ARTS_USER_ERROR_IF(group_index < 0, "Cannot recognize group: ", group)

        WsvRecord x(name, desc ? desc : "Created by pybind11 API", group_index);

        Index out = w.add_wsv_inplace(x);
        w.push(out, workspace_memory_handler.allocate(group_index));
      },
      py::arg("group"),
      py::arg("name"),
      py::arg("desc") = nullptr,
      py::doc(
          R"--(
Creates a named variable of the given group, initializing it on the workspace

The variable can be accessed with its property name upon completion

Example:
ws.create_variable("String", "my_string")

Now,
ws.my_string

can be accessed and will give an empty string object.  This you can
ws.my_string = "hi there"

and now doing
print(ws.my_string)

will output the expected greeting
)--"));

  ws.def("__delattr__", [](Workspace& w, const char* name) {
    auto varpos = std::find_if(Workspace::WsvMap.begin(),
                               Workspace::WsvMap.end(),
                               [name](auto& b) { return b.first == name; });
    ARTS_USER_ERROR_IF(varpos == Workspace::WsvMap.end(),
                       "No workspace variable: ",
                       name,
                       "\n\nCustom workspace variables have to "
                       "be created using create_variable or by explicit set")

    if (w.is_initialized(varpos->second)) {
      w.pop_free(varpos->second);
    }
  });

  ws.def("__getattr__", [](Workspace& w, const char* name) {
    auto varpos = std::find_if(Workspace::WsvMap.begin(),
                               Workspace::WsvMap.end(),
                               [name](auto& b) { return b.first == name; });
    ARTS_USER_ERROR_IF(varpos == Workspace::WsvMap.end(),
                       "No workspace variable: ",
                       name,
                       "\n\nCustom workspace variables have to "
                       "be created using create_variable or by explicit set")

    return WorkspaceVariable{w, varpos->second};
  });

  ws.def(
      "swap",
      [](Workspace* w_, Workspace* w_2) { std::swap(w_, w_2); },
      py::return_value_policy::reference_internal);

  py::class_<WorkspaceVariable>(m, "WorkspaceVariable")
      .def_property(
          "value",
          py::cpp_function([](const WorkspaceVariable& wsv)
                               -> WorkspaceVariablesVariant { return wsv; },
                           py::return_value_policy::reference_internal),
          [](WorkspaceVariable& wsv, WorkspaceVariablesVariant x) { wsv = x; },
          "Returns a proper Arts type")
          .def_property_readonly("name", &WorkspaceVariable::name)
      .def("__str__",
           [](WorkspaceVariable& wsv) {
             return var_string(
                 "Arts Workspace Variable ",
                 wsv.ws.wsv_data[wsv.pos].Name(),
                 " of type ",
                 global_data::wsv_group_names[Workspace::wsv_data[wsv.pos]
                                                  .Group()]);
           })
      .def("__repr__", [](WorkspaceVariable& wsv) {
        return var_string(
            "Arts Workspace Variable ",
            wsv.ws.wsv_data[wsv.pos].Name(),
            " of type ",
            global_data::wsv_group_names[Workspace::wsv_data[wsv.pos].Group()]);
      });

  // Should be last as it contains several implicit conversions
  py_auto_workspace(ws);
}
}  // namespace Python