#include <auto_md.h>
#include <global_data.h>
#include <py_auto_interface.h>
#include <xml_io.h>

#include <algorithm>
#include <filesystem>

#include "py_macros.h"
#include <parameters.h>

extern Parameters parameters;

namespace Python {
Index create_workspace_gin_default_internal(Workspace& ws, String& key);

Agenda *parse_agenda(const char *filename, const Verbosity& verbosity) {
  Agenda *a = new Agenda;
  ArtsParser parser = ArtsParser(*a, filename, verbosity);

  parser.parse_tasklist();
  a->set_name(filename);
  a->set_main_agenda();
  return a;
}

struct AgendaMethodVariable {
  String name{};
  Index group{-1};
  Index ws_pos{-1};
  Index g_pos{-1};
  bool input{false};
  bool any{false};
};

Array<AgendaMethodVariable> sorted_mdrecord(
    Workspace& ws, decltype(global_data::md_data_raw.begin()) ptr) {
  static std::map<String, Index> gin_defaults;

  const Index any_pos = global_data::WsvGroupMap.at("Any");

  Array<AgendaMethodVariable> var_order;

  for (auto& output : ptr->Out()) {
    auto& var = var_order.emplace_back();
    var.name = ws.wsv_data[output].Name();
    var.group = ws.wsv_data[output].Group();
    var.ws_pos = output;
    var.input = false;
  }

  for (Index i = 0; i < ptr->GOut().nelem(); i++) {
    auto& var = var_order.emplace_back();
    var.name = ptr->GOut()[i];
    var.group = ptr->GOutType()[i];
    var.g_pos = i;
    var.input = false;
    var.any = var.group == any_pos;
  }

  for (auto& input : ptr->In()) {
    // Skip inout
    if (std::any_of(ptr->Out().begin(),
                    ptr->Out().end(),
                    [input](auto& output) { return input == output; }))
      continue;

    auto& var = var_order.emplace_back();
    var.name = ws.wsv_data[input].Name();
    var.group = ws.wsv_data[input].Group();
    var.ws_pos = input;
    var.input = true;
  }

  for (Index i = 0; i < ptr->GIn().nelem(); i++) {
    auto& var = var_order.emplace_back();
    var.name = ptr->GIn()[i];
    var.group = ptr->GInType()[i];
    var.g_pos = i;
    var.input = true;
    var.any = var.group == any_pos;

    // Create defaults (on the workspace)
    if (ptr->GInDefault()[i] not_eq NODEF) {
      String gin_key =
          "GeneratedInputDefault_" + ptr->Name() + "_" + var.name + "_";

      auto gin_def_ptr = gin_defaults.find(gin_key);
      if (gin_def_ptr == gin_defaults.end()) {
        gin_defaults[gin_key] =
            create_workspace_gin_default_internal(ws, gin_key);
        gin_def_ptr = gin_defaults.find(gin_key);
      }

      var.ws_pos = gin_def_ptr->second;
    }
  }

  return var_order;
}

std::pair<ArrayOfIndex, ArrayOfIndex> split_io(
    Array<AgendaMethodVariable>& var_order) {
  ArrayOfIndex in, out;
  for (auto& var : var_order) {
    if (var.input)
      in.push_back(var.ws_pos);
    else
      out.push_back(var.ws_pos);
  }
  return {in, out};
}

void py_agenda(py::module_& m) {
  py::class_<CallbackFunction>(m, "CallbackFunction")
      .def(py::init<>())
      .def(py::init<std::function<void(Workspace&)>>())
      .def("__call__", [](CallbackFunction& f, Workspace& ws) { f(ws); });
  py::implicitly_convertible<std::function<void(Workspace&)>,
                             CallbackFunction>();

  py::class_<MdRecord>(m, "MdRecord")
      .def(py::init<>())
      .def_property_readonly("name", &MdRecord::Name)
      .def_property_readonly("out", &MdRecord::Out)
      .def_property_readonly("gout", &MdRecord::GOut)
      .def_property_readonly("gout_type", &MdRecord::GOutType)
      .def_property_readonly("in", &MdRecord::In)
      .def_property_readonly("gin", &MdRecord::GIn)
      .def_property_readonly("gin_type", &MdRecord::GInType);

  py::class_<Array<MdRecord>>(m, "ArrayOfMdRecord")
      .PythonInterfaceArrayDefault(MdRecord);

  py::class_<Agenda>(m, "Agenda")
      .def(py::init<>())
      .def(py::init([](Workspace& w, std::filesystem::path path) {
        ARTS_USER_ERROR_IF(not std::filesystem::is_regular_file(path), "There is no regular file at: ", path.c_str())

        std::filesystem::path path_copy = path;
        parameters.includepath.push_back(path_copy.remove_filename().c_str());
        parameters.includepath.erase(std::unique(parameters.includepath.begin(), parameters.includepath.end()), parameters.includepath.end());

        Agenda *a = parse_agenda(path.c_str(), *reinterpret_cast<Verbosity *>(w[w.WsvMap.at("verbosity")]));
        w.initialize();
        return a;
      }))
      .def(
          "add_workspace_method",
          [](Agenda& a,
             Workspace& ws,
             const char* name,
             py::args args,
             py::kwargs kwargs) {
            const Index nargs = args.size();

            // The MdRecord
            auto ptr = std::find_if(
                global_data::md_data_raw.begin(),
                global_data::md_data_raw.end(),
                [name](auto& method) { return name == method.Name(); });
            ARTS_USER_ERROR_IF(ptr == global_data::md_data_raw.end(),
                               "No \"",
                               name,
                               "\" method found")

            // The method input and output lists in one long list
            Array<AgendaMethodVariable> var_order = sorted_mdrecord(ws, ptr);

            // Generic checker
            const bool generic_function =
                std::any_of(var_order.begin(), var_order.end(), [](auto& v) {
                  return v.any;
                });

            // Point at the first actual function (or if non-generic, the actual function)
            ptr = std::find_if(
                global_data::md_data.begin(),
                global_data::md_data.end(),
                [name](auto& method) { return name == method.Name(); });

            // Set the actual input and output of all the variables (and create anonymous variables)
            for (Index i = 0; i < var_order.nelem(); i++) {
              auto& var = var_order[i];
              try {
                if (i < nargs) {
                  var.ws_pos = WorkspaceVariable(ws, var.group, args[i]).pos;
                } else if (auto key = var.name.c_str(); kwargs.contains(key)) {
                  var.ws_pos =
                      WorkspaceVariable(ws, var.group, kwargs[key]).pos;
                }
              } catch (std::exception& e) {
                ARTS_USER_ERROR("\nCannot add workspace method \"",
                                name,
                                "\" with the following reason for variable \"",
                                var.name,
                                "\":\n\n",
                                e.what());
              }
            }

            // Adjust pointer to the correct one given the current args and kwargs if we are generic
            if (generic_function) {
              // Ensure we have GIN and GOUT for all of the generic types
              for (Index i = 0; i < var_order.nelem(); i++) {
                auto& var = var_order[i];
                ARTS_USER_ERROR_IF(var.ws_pos < 0,
                                   "Must have variable ",
                                   var.name,
                                   " for agenda method ",
                                   name)
                var.group = ws.wsv_data[var.ws_pos].Group();
              }

              // Loop until all generic GIN and GOUT fits a known method
              while (ptr not_eq global_data::md_data.end() and
                     ptr->Name() == name and
                     std::any_of(
                         var_order.begin(), var_order.end(), [ptr](auto& var) {
                           return var.any and
                                  (var.input ? ptr->GInType()[var.g_pos] not_eq
                                                   var.group
                                             : ptr->GOutType()[var.g_pos] not_eq
                                                   var.group);
                         })) {
                ptr++;
              }

              ARTS_USER_ERROR_IF(
                  ptr == global_data::md_data.end() or ptr->Name() not_eq name,
                  "Cannot find a matching generic function signature for call to: ",
                  name,
                  "\n"
                  "Note that generic functions can only operate on true Workspace groups, "
                  "either via a workspace variable or via python classes representing the type")
            }

            // Ensure we have all input/output
            for (auto& var : var_order) {
              ARTS_USER_ERROR_IF(var.ws_pos < 0,
                                 "Must set ",
                                 var.name,
                                 " in agenda method ",
                                 name)
            }

            // Set the method record
            const auto [in, out] = split_io(var_order);
            a.push_back(
                MRecord(std::distance(global_data::md_data.begin(), ptr),
                        out,
                        in,
                        {},
                        {}));
          },
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
          [](Agenda& a, Workspace& ws, CallbackFunction x) mutable {
            static std::size_t counter = 0;

            const Index group_index = get_wsv_group_id("CallbackFunction");
            ARTS_USER_ERROR_IF(group_index < 0,
                               "Cannot recognize CallbackFunction")

            const String name = var_string("_callback_function_nr_", counter++);
            Index in = ws.add_wsv_inplace(WsvRecord(
                name.c_str(), "Callback created by pybind11 API", group_index));
            ws.push(in, new CallbackFunction{x});

            auto method_ptr =
                global_data::MdMap.find("CallbackFunctionExecute");
            ARTS_USER_ERROR_IF(method_ptr == global_data::MdMap.end(),
                               "Cannot find CallbackFunctionExecute")

            a.push_back(MRecord(method_ptr->second, {}, {in}, {}, {}));
          })
      .def("append_agenda_methods", [](Agenda& self, Agenda& other) {
        for(auto& method: other.Methods()) self.push_back(method);
      })
      .def("execute", &Agenda::execute)
      .def("check", [](Agenda& a, Workspace& w) { a.check(w, *reinterpret_cast<Verbosity *>(w[w.WsvMap.at("verbosity")])); })
      .def_property("name", &Agenda::name, &Agenda::set_name)
      .def("__repr__",
           [](Agenda& a) {
             std::string out = var_string("Arts Agenda ", a.name(), ":\n");
             std::ostringstream os;
             a.print(os, "   ");
             out += os.str();
             return out;
           })
      .def("__str__", [](Agenda& a) {
        std::string out = var_string("Arts Agenda ", a.name(), ":\n");
        std::ostringstream os;
        a.print(os, "   ");
        out += os.str();
        return out;
      });

  py::class_<ArrayOfAgenda>(m, "ArrayOfAgenda");
}
}  // namespace Python