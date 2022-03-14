#include <auto_md.h>
#include <global_data.h>
#include <parameters.h>
#include <py_auto_interface.h>
#include <pybind11/cast.h>
#include <xml_io.h>

#include <algorithm>
#include <filesystem>
#include <iterator>

#include "agenda_class.h"
#include "array.h"
#include "debug.h"
#include "py_macros.h"
#include "tokval.h"

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

Agenda* parse_agenda(const char* filename, const Verbosity& verbosity) {
  auto* a = new Agenda;
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
      var.ws_pos = create_workspace_gin_default_internal(
          ws, var_string("::", ptr->Name(), "::", var.name));
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

MRecord simple_set_method(WorkspaceVariable val) {
  using namespace global_data;

  // Find group and is it acceptable
  const String group = wsv_group_names[val.group()];
  ARTS_USER_ERROR_IF(group == "ArrayOfAgenda",
                     "Cannot support setting ArrayOfAgenda")

  // Set method set values
  TokVal t;
  Agenda a;
  if (group == "Agenda") {
    a = Agenda(val);
  } else {
    t = TokVal(val);
  }

  // If this is a new variable, we need to find its set-method and transfer
  // data ownership to the MRecord and remove the value from the workspace
  const ArrayOfIndex output = ArrayOfIndex{val.pos};
  const Index m_id = MdMap.find(group + "Set")->second;
  val.pop_workspace_level();

  return MRecord(m_id, output, {}, t, a);
}

MRecord simple_delete_method(WorkspaceVariable val) {
  using namespace global_data;

  // Find group and is it acceptable
  auto ptr = std::find_if(md_data.begin(),
                          md_data.end(),
                          [delete_method = String("Delete")](auto& m) {
                            return m.Name() == delete_method;
                          });
  while (ptr not_eq md_data.end() and ptr->In().nelem() == 1 and
         ptr->In().front() not_eq val.group())
    ptr++;
  ARTS_ASSERT(ptr not_eq md_data.end(),
              "Did not find a Delete method to the end of the list")
  ARTS_ASSERT(ptr->Name() == "Delete",
              "Expects to find a Delete method but encounters ",
              ptr->Name())

  return MRecord(std::distance(md_data.begin(), ptr), {}, {val.pos}, {}, {});
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
      .def(py::init(
          [](Workspace& ws, const std::function<py::object(Workspace&)>& f) {
            return py::cast<Agenda>(f(ws));
          }))
      .PythonInterfaceWorkspaceVariableConversion(Agenda)
      .def(py::init([](Workspace& w, const std::filesystem::path& path) {
        Agenda* a = parse_agenda(
            correct_include_path(path).c_str(),
            *reinterpret_cast<Verbosity*>(w[w.WsvMap.at("verbosity")]));
        w.initialize();
        return a;
      }))
      .PythonInterfaceFileIO(Agenda)
      .def_property_readonly("main", &Agenda::is_main_agenda)
      .def_property("name", &Agenda::name, &Agenda::set_name)
      .def_property_readonly("output2push", &Agenda::get_output2push)
      .def_property_readonly("output2dup", &Agenda::get_output2dup)
      .def("set_outputs_to_push_and_dup",
           [](Agenda& a, Workspace& w) {
             a.set_outputs_to_push_and_dup(
                 *reinterpret_cast<Verbosity*>(w[w.WsvMap.at("verbosity")]));
           })
      .def(
          "add_workspace_method",
          [](Agenda& a,
             Workspace& ws,
             const char* name,
             const py::args& args,
             const py::kwargs& kwargs) {
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
                } else if (const auto* key = var.name.c_str();
                           kwargs.contains(key)) {
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

            // Set the actual values using automatic set methods
            for (auto& var : var_order) {
              WorkspaceVariable val(ws, var.ws_pos);
              if (val.name().find("::anon::") == 0) {
                a.push_back(simple_set_method(val));
              }
            }

            // Set the method record
            const auto [in, out] = split_io(var_order);
            const auto m_id = std::distance(global_data::md_data.begin(), ptr);
            if (ptr->SetMethod() and in.nelem()) {
              ARTS_USER_ERROR_IF (ws.is_initialized(in.back()), "Can only set values, use Copy to copy workspace variables")

              // The previous value will be the GIN set value as the last value
              // is also the last in the push_back above
              ARTS_ASSERT(a.Methods().nelem())
              a.push_back(MRecord(m_id, out, in, a.Methods().back().SetValue(), a.Methods().back().Tasks()));
            } else {
              a.push_back(MRecord(m_id, out, in, {}, {}));
            }

            // Clean up the value
            for (auto& var : var_order) {
              WorkspaceVariable val(ws, var.ws_pos);
              if (val.name().find("::anon::") == 0) {
                a.push_back(simple_delete_method(val));
              }
            }
          },
          py::arg("ws"),
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
          [](Agenda& a, Workspace& ws, const CallbackFunction& f) mutable {
            const Index group_index = get_wsv_group_id("CallbackFunction");
            ARTS_USER_ERROR_IF(group_index < 0,
                               "Cannot recognize CallbackFunction")

            static std::size_t counter = 0;
            const String name = var_string("::callback::", counter++);

            Index in = ws.add_wsv_inplace(WsvRecord(
                name.c_str(), "Callback created by pybind11 API", group_index));
            ws.push(in, new CallbackFunction{f});

            auto method_ptr =
                global_data::MdMap.find("CallbackFunctionExecute");
            ARTS_USER_ERROR_IF(method_ptr == global_data::MdMap.end(),
                               "Cannot find CallbackFunctionExecute")

            WorkspaceVariable val(ws, in);
            a.push_back(simple_set_method(val));
            a.push_back(MRecord(method_ptr->second, {}, {in}, {}, {}));
            a.push_back(simple_delete_method(val));
          },
          py::doc(R"--(
Adds a callback method to the Agenda

Parameters
----------
    ws : Workspace
        The workspace for which the Agenda works
    f : CallbackFunction
        A method that takes ws and only ws as input
)--"),
          py::arg("ws"),
          py::arg("f"))
      .def(
          "append_agenda_methods",
          [](Agenda& self, const Agenda& other) {
            const Index n = other.Methods().nelem();  // if other==self
            for (Index i=0 ; i<n; i++) self.push_back(other.Methods()[i]);
          },
          py::doc(R"--(
Appends the input agenda's methods to this agenda's method list

Parameters
----------
other : Agenda
    The other agenda

Warning
-------
Both agendas must be defined on the same workspace)--"), py::arg("other"))
      .def(
          "execute",
          [](Agenda& a, Workspace& ws) {
            a.set_main_agenda();
            a.execute(ws);
          },
          "Executes the agenda as if it was the main agenda")
      .def(
          "check",
          [](Agenda& a, Workspace& w) {
            a.check(w,
                    *reinterpret_cast<Verbosity*>(w[w.WsvMap.at("verbosity")]));
          },
          py::doc("Checks if the agenda works"))
      .def_property("name", &Agenda::name, &Agenda::set_name)
      .def("__repr__",
           [](Agenda& a) { return var_string("Agenda ", a.name()); })
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
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfAgenda)
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
      .def("__repr__", [](ArrayOfAgenda&) { return "ArrayOfAgenda"; })
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
      .def("pop", [](ArrayOfAgenda& aa) {
        Agenda a = std::move(aa.back());
        aa.pop_back();
        return a;
      }, "Pops and returns the last item")
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
