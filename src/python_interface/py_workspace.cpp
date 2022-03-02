
#include <depr.h>
#include <global_data.h>
#include <parameters.h>
#include <py_auto_interface.h>
#include <workspace_ng.h>

#include <cstdlib>

#include "debug.h"
#include "py_macros.h"

extern Parameters parameters;

void parse_path_from_environment(String envvar, ArrayOfString& paths);
namespace Python {
std::filesystem::path correct_include_path(
    const std::filesystem::path& path_copy);
Agenda* parse_agenda(const char* filename, const Verbosity& verbosity);
void py_auto_workspace(py::class_<Workspace>&);

void py_workspace(py::module_& m, py::class_<Workspace>& ws) {
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
           Workspace w{};
           w.initialize();
           w.push(w.WsvMap.at("verbosity"),
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
             const char* name,
             std::optional<const char*> desc) {
            using global_data::workspace_memory_handler;

            const Index group_index = get_wsv_group_id(group);
            ARTS_USER_ERROR_IF(
                group_index < 0, "Cannot recognize group: ", group)

            auto varptr = w.WsvMap.find(name);
            ARTS_USER_ERROR_IF(
                varptr not_eq w.WsvMap.end() and
                    w.wsv_data[varptr->second].Group() not_eq group_index,
                "Already exist of different group: ",
                name)

            w.push(w.add_wsv_inplace(
                       WsvRecord(name,
                                 desc.has_value() ? *desc : "No description",
                                 group_index)),
                   workspace_memory_handler.allocate(group_index));
          },
          py::arg("group").none(false),
          py::arg("name").none(false),
          py::arg("desc") = std::nullopt)
      .def(
          "__delattr__",
          [](Workspace& w, const char* name) {
            auto varpos = w.WsvMap.find(name);
            if (varpos == w.WsvMap.end())
              throw pybind11::attribute_error(var_string('\'', name, '\''));

            if (w.is_initialized(varpos->second)) {
              w.pop_free(varpos->second);
            }
          },
          py::arg("name").none(false));

  ws.def(
      "__getattr__",
      [](Workspace& w, const char* name) {
        auto varpos = w.WsvMap.find(name);
        if (varpos == w.WsvMap.end())
          throw pybind11::attribute_error(var_string('\'', name, '\''));

        return WorkspaceVariable{w, varpos->second};
      },
      py::doc("A custom workspace variable"));

  ws.def(
      "swap",
      [](Workspace& w_1, Workspace& w_2) { w_1.swap(w_2); },
      py::arg("ws").noconvert().none(false));

  ws.def("__repr__", [](Workspace&) { return "Workspace"; });

  ws.def("__str__", [](Workspace& w) {
    Index c = 0;
    for (Index i = 0; i < w.nelem(); i++) c += w.is_initialized(i);
    return var_string("Workspace with ", c, " initialized variables");
  });

  py::class_<WorkspaceVariable>(m, "WorkspaceVariable")
      .def_property(
          "value",
          py::cpp_function([](const WorkspaceVariable& wsv)
                               -> WorkspaceVariablesVariant { return wsv; },
                           py::return_value_policy::reference_internal),
          [](WorkspaceVariable& wsv, WorkspaceVariablesVariant x) { wsv = x; },
          "Returns a proper Arts type")
      .def_property_readonly("name", &WorkspaceVariable::name)
      .def_property_readonly("desc", [](WorkspaceVariable& wv) {return wv.ws.wsv_data[wv.pos].Description();})
      .def_property_readonly("init", &WorkspaceVariable::is_initialized)
      .def("initialize_if_not", &WorkspaceVariable::initialize_if_not)
      .def(
          "__str__",
          [](WorkspaceVariable& wsv) {
            return var_string(
                "Workspace ",
                global_data::wsv_group_names[wsv.ws.wsv_data[wsv.pos].Group()]);
          })
      .def("__repr__", [](WorkspaceVariable& wsv) {
        return var_string(
            "Workspace ",
            global_data::wsv_group_names[wsv.ws.wsv_data[wsv.pos].Group()]);
      });

  py::class_<WsvRecord>(m, "WsvRecord")
      .def_property_readonly("name", &WsvRecord::Name)
      .def_property_readonly("description", &WsvRecord::Description)
      .def_property_readonly("group", &WsvRecord::Group)
      .def_property_readonly("implicit", &WsvRecord::Implicit)
      .def("__str__", [](WsvRecord& mr) { return var_string(mr); })
      .def("__repr__", [](WsvRecord& mr) { return var_string(mr); });

  py::class_<Array<WsvRecord>>(m, "ArrayOfWsvRecord")
      .PythonInterfaceArrayDefault(WsvRecord)
      .PythonInterfaceBasicRepresentation(Array<WsvRecord>);

  // Should be last as it contains several implicit conversions
  py_auto_workspace(ws);
}
}  // namespace Python