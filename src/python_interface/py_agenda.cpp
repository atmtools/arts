#include <auto_md.h>
#include <global_data.h>
#include <parameters.h>
#include <py_auto_interface.h>
#include <xml_io.h>

#include <algorithm>
#include <filesystem>

#include "py_macros.h"

extern Parameters parameters;

namespace Python {
Index create_workspace_gin_default_internal(Workspace& ws, const String& key);

std::filesystem::path correct_include_path(std::filesystem::path path) {
  const std::filesystem::path path_copy=path;
  for (auto& prefix : parameters.includepath) {
    if (std::filesystem::is_regular_file(
            path = std::filesystem::path(prefix.c_str()) / path_copy))
      break;
  }

  ARTS_USER_ERROR_IF(not std::filesystem::is_regular_file(path),
                     "There is no regular file at ",
                     path_copy,
                     '\n',
                     "Looked in the following search path ",
                     '[',
                     parameters.includepath,
                     ']')

  return path;
}

Agenda* parse_agenda(const char* filename, const Verbosity& verbosity) {
  Agenda* a = new Agenda;
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
      const String gin_key = var_string("::", ptr->Name(), "::", var.name);

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

String error_msg(const String& name,
                 Workspace& ws,
                 const Array<AgendaMethodVariable>& var_order) {
  std::ostringstream os;
  os << name << "(";
  for (auto& var : var_order) {
    os << var.name << " : ";
    if (var.ws_pos >= 0)
      os << global_data::wsv_group_names[WorkspaceVariable(ws, var.ws_pos)
                                             .group()];
    else
      os << global_data::wsv_group_names[var.group];
    if (&var not_eq &var_order.back()) os << ", ";
  }
  os << ")";
  return os.str();
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
      .def_property_readonly("outs", &MdRecord::Out)
      .def_property_readonly("g_out", &MdRecord::GOut)
      .def_property_readonly("g_out_types", &MdRecord::GOutType)
      .def_property_readonly("ins", &MdRecord::In)
      .def_property_readonly("g_in", &MdRecord::GIn)
      .def_property_readonly("g_in_types", &MdRecord::GInType)
      .def_property_readonly("g_in_default", &MdRecord::GInDefault)
      .def("__str__", [](MdRecord& mr) { return var_string('"', mr, '"'); })
      .def("__repr__", [](MdRecord& mr) { return var_string('"', mr, '"'); });

  py::class_<Array<MdRecord>>(m, "ArrayOfMdRecord")
      .PythonInterfaceArrayDefault(MdRecord)
      .PythonInterfaceBasicRepresentation(Array<MdRecord>);

  py::class_<AgRecord>(m, "AgRecord")
      .def(py::init<>())
      .def_property_readonly("name", &AgRecord::Name)
      .def_property_readonly("description", &AgRecord::Description)
      .def_property_readonly("outs", &AgRecord::Out)
      .def_property_readonly("ins", &AgRecord::In)
      .def("__str__", [](AgRecord& ar) { return var_string(ar); })
      .def("__repr__", [](AgRecord& ar) { return var_string(ar); });

  py::class_<Array<AgRecord>>(m, "ArrayOfAgRecord")
      .PythonInterfaceArrayDefault(AgRecord)
      .PythonInterfaceBasicRepresentation(Array<AgRecord>);

  py::class_<Agenda>(m, "Agenda")
      .def(py::init<>())
      .def(py::init([](Workspace& w, std::filesystem::path path) {
        Agenda* a = parse_agenda(
            correct_include_path(path).c_str(),
            *reinterpret_cast<Verbosity*>(w[w.WsvMap.at("verbosity")]));
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

            try_to_set_workspace_position:
              try {
                if (i < nargs) {
                  var.ws_pos =
                      WorkspaceVariable(ws, var.group, args[i], var.input).pos;
                } else if (auto key = var.name.c_str(); kwargs.contains(key)) {
                  var.ws_pos =
                      WorkspaceVariable(ws, var.group, kwargs[key], var.input)
                          .pos;
                }
              } catch (std::exception& e) {
                // We should be able to switch to the right method for simple_generic_function
                if (var.any and var.input) {
                increment_ptr:
                  ptr++;

                  ARTS_USER_ERROR_IF(
                      ptr == global_data::md_data.end() or
                          ptr->Name() not_eq name,
                      "Cannot find matching workspace method:\n",
                      error_msg(name, ws, var_order),
                      "\n"
                      "The generic variable \"",
                      var.name,
                      "\" cannot be converted to a type that matches any available workspace method signatures (note that the signature above is therefore potentially extra weird)")

                  // We have a new variable order but previous variables must be respected
                  auto new_var_order = sorted_mdrecord(ws, ptr);
                  for (Index j = 0; j < i; j++) {
                    if (var_order[j].ws_pos >= 0 and
                        new_var_order[j].group not_eq
                            WorkspaceVariable(ws, var_order[j].ws_pos)
                                .group()) {
                      goto increment_ptr;
                    }
                  }

                  // Try with a new type
                  var.group = new_var_order[i].group;
                  goto try_to_set_workspace_position;
                } else {
                  ARTS_USER_ERROR("\nCannot add workspace method:\n",
                                  error_msg(name, ws, var_order),
                                  "\n"
                                  "with the following reason for variable \"",
                                  var.name,
                                  "\":\n\n",
                                  e.what());
                }
              }
            }

            // Ensure we have all input/output
            for (auto& var : var_order) {
              ARTS_USER_ERROR_IF(var.ws_pos < 0,
                                 "Must set ",
                                 var.name,
                                 " in agenda method:\n",
                                 error_msg(name, ws, var_order))
            }

            // Adjust pointer to the correct one given the current args and kwargs if we are generic
            if (generic_function) {
              // Ensure we have GIN and GOUT for all of the generic types
              for (Index i = 0; i < var_order.nelem(); i++) {
                auto& var = var_order[i];
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
                  "Cannot find a matching generic function signature for call to:\n",
                  error_msg(name, ws, var_order),
                  "\n"
                  "Please check documentation of the function for correct signatures\n\n"
                  "If you think we should have deduced a correct signature from your\n"
                  "input, try first casting your data into a Arts type manually or put\n"
                  "it directly onto the workspace")
            } else {
              // Ensure that all input match the final method
              ARTS_USER_ERROR_IF(
                  std::any_of(
                      var_order.begin(),
                      var_order.end(),
                      [&](auto& var) {
                        return WorkspaceVariable(ws, var.ws_pos).group() not_eq
                               var.group;
                      }),
                  "Mismatching signature.  Got:\n",
                  error_msg(name, ws, var_order),
                  "\n"
                  "But this does not match a signature in Arts")
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
            const Index group_index = get_wsv_group_id("CallbackFunction");
            ARTS_USER_ERROR_IF(group_index < 0,
                               "Cannot recognize CallbackFunction")

            static std::size_t counter = 0;
            const String name = var_string("::callback::", counter++);

            Index in = ws.add_wsv_inplace(WsvRecord(
                name.c_str(), "Callback created by pybind11 API", group_index));
            ws.push(in, new CallbackFunction{x});

            auto method_ptr =
                global_data::MdMap.find("CallbackFunctionExecute");
            ARTS_USER_ERROR_IF(method_ptr == global_data::MdMap.end(),
                               "Cannot find CallbackFunctionExecute")

            a.push_back(MRecord(method_ptr->second, {}, {in}, {}, {}));
          })
      .def("append_agenda_methods",
           [](Agenda& self, Agenda& other) {
             for (auto& method : other.Methods()) self.push_back(method);
           })
      .def("execute", &Agenda::execute)
      .def(
          "check",
          [](Agenda& a, Workspace& w) {
            a.check(w,
                    *reinterpret_cast<Verbosity*>(w[w.WsvMap.at("verbosity")]));
          },
          py::doc("Checks if the agenda works"))
      .def_property("name", &Agenda::name, &Agenda::set_name)
      .def("__repr__",
           [](Agenda& a) { return var_string("Arts Agenda ", a.name()); })
      .def("__str__", [](Agenda& a) {
        std::string out =
            var_string("Arts Agenda ", a.name(), " with methods:\n");
        std::ostringstream os;
        a.print(os, "   ");
        out += os.str();
        return out;
      });

  py::class_<ArrayOfAgenda>(m, "ArrayOfAgenda")
      .def(py::init<>())
      .def(py::init<Index>())
      .def(py::init<Index, Agenda>())
      .def(py::init([](std::vector<Agenda> va) {
        for (auto& a : va) {
          ARTS_USER_ERROR_IF(
              a.name() not_eq va.front().name(),
              "An ArrayOfAgenda must only consist of agendas with the same name\n"
              "You have input a list of agendas that contains disimilar names.\n"
              "\nThe first item is named: \"",
              va.front().name(),
              '"',
              '\n',
              "A later item in the list is names: \"",
              a.name(),
              '"',
              '\n')
        }
        return va;
      }))
      .def("__repr__",
           [](ArrayOfAgenda& aa) {
             std::string out = "[\n";
             for (auto& a : aa) {
               out += var_string("Arts Agenda ", a.name(), ":\n");
               std::ostringstream os;
               a.print(os, "   ");
               out += os.str();
               out += "\n";
             }
             out += "]\n";
             return out;
           })
      .def("__str__",
           [](ArrayOfAgenda& aa) {
             std::string out = "[\n";
             for (auto& a : aa) {
               out += var_string("Arts Agenda ", a.name(), ":\n");
               std::ostringstream os;
               a.print(os, "   ");
               out += os.str();
               out += "\n";
             }
             out += "]\n";
             return out;
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
            if (y.name() not_eq x.front().name()) y.set_name(x.front().name());
            x[i] = std::move(y);
          },
          py::return_value_policy::reference_internal)
      .def("append",
           [](ArrayOfAgenda& x, Agenda y) {
             if (x.nelem() and y.name() not_eq x.front().name())
               y.set_name(x.front().name());
             x.emplace_back(std::move(y));
           })
      .def(
          "check",
          [](ArrayOfAgenda& aa, Workspace& w) {
            for (auto& a : aa)
              a.check(
                  w,
                  *reinterpret_cast<Verbosity*>(w[w.WsvMap.at("verbosity")]));
          },
          py::doc("Checks if the agenda works"))
      .def_property(
          "name",
          [](ArrayOfAgenda& a) -> String {
            if (a.nelem() == 0) return "";
            return a.front().name();
          },
          [](ArrayOfAgenda& aa, const String& name) {
            for (auto& a : aa) a.set_name(name);
          });
  py::implicitly_convertible<std::vector<Agenda>, ArrayOfAgenda>();
}
}  // namespace Python
