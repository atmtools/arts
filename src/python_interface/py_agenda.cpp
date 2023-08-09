#include <algorithm>
#include <global_data.h>
#include <parameters.h>
#include <parser.h>
#include <pybind11/attr.h>
#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl/filesystem.h>
#include <memory>
#include <type_traits>

#include "agenda_class.h"
#include "debug.h"
#include "py_auto_interface.h"
#include "py_macros.h"
#include "python_interface.h"
#include "tokval_variant.h"
#include "wsv_aux.h"

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

std::unique_ptr<Agenda> parse_agenda(Workspace& ws, const char* filename, const Verbosity& verbosity) {
  std::unique_ptr<Agenda> a = std::make_unique<Agenda>(ws);
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
    var.name = ws.wsv_data_ptr->operator[](output).Name();
    var.group = ws.wsv_data_ptr->operator[](output).Group();
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
    var.name = ws.wsv_data_ptr->operator[](input).Name();
    var.group = ws.wsv_data_ptr->operator[](input).Group();
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
      os << global_data::wsv_groups[WorkspaceVariable(ws, var.ws_pos)
                                             .group()];
    else
      os << global_data::wsv_groups[var.group];
    if (&var not_eq &var_order.back()) os << ", ";
  }
  os << ")";
  return os.str();
}

MRecord simple_set_method(WorkspaceVariable val) {
  using namespace global_data;

  // Find group and is it acceptable
  const auto& group = wsv_groups[val.group()];
  ARTS_USER_ERROR_IF(group == "ArrayOfAgenda",
                     "Cannot support setting ArrayOfAgenda")

  // Set method set values
  TokVal t;
  Agenda a{val.ws};
  if (group == "Agenda") {
    a = Agenda(val);
  } else {
    t = TokVal(val);
  }

  // If this is a new variable, we need to find its set-method and transfer
  // data ownership to the MRecord and remove the value from the workspace
  const ArrayOfIndex output = ArrayOfIndex{val.pos};
  const Index m_id = MdMap.find(group.name + "Set")->second;
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

  return MRecord(
      std::distance(md_data.begin(), ptr),
      {},
      {val.pos},
      {},
      Agenda{val.ws});
}

void py_agenda(py::module_& m) {
  py::class_<CallbackFunction>(m, "CallbackFunction")
      .def(py::init([]() { return std::make_unique<CallbackFunction>(); }), py::doc("Initialize as empty call"))
      .PythonInterfaceCopyValue(CallbackFunction)
      .def(py::init([](const std::function<void(Workspace&)>& f) {
        return std::make_unique<CallbackFunction>(f);
      }), py::doc("Initialize from python callable"))
      .def(
          "__call__",
          [](CallbackFunction& f, Workspace& ws) { f(ws); })
      .PythonInterfaceWorkspaceDocumentationExtra(CallbackFunction, R"(

An instance of this object can be called taking only an
instance of the workspace as input.

Example
-------
>>> import pyarts
>>> def print_ws(ws):
>>>     print(ws.iy_unit.value)
>>> ws = pyarts.workspace.Workspace()
>>> cb = pyarts.arts.CallbackFunction(print_ws)
>>> cb(ws)
1
)");
  py::implicitly_convertible<std::function<void(Workspace&)>,
                             CallbackFunction>();

  py::class_<TokVal>(m, "TokVal")
      .def(py::init([](const py::int_& x) { return TokVal{x.cast<Index>()}; }), "Holds :class:`~pyarts.arts.Index` from :class:`int`")
      .def(py::init([](const py::float_& x) { return TokVal{x.cast<Numeric>()}; }), "Holds :class:`~pyarts.arts.Numeric` from :class:`float`")
      .def(py::init([](const WorkspaceVariablesVariant& x) -> TokVal {
        return std::visit(
            [](auto& v) -> TokVal {
              if constexpr (std::is_same_v<Index_*,
                                           std::remove_cvref_t<decltype(v)>> or
                            std::is_same_v<Numeric_*,
                                           std::remove_cvref_t<decltype(v)>>)
                return v->val;
              else if constexpr (
                  std::is_same_v<py::object*, std::remove_cvref_t<decltype(v)>>)
                ARTS_USER_ERROR(
                    "Must be a non-agenda ARTS type (try specifying the type if the input is not an Agenda type)")
              else
                return *v;
            },
            x);
      }), "Holds custom ARTS Workspace group")
      .def_property(
          "value",
          [](const TokVal& x) {
            return std::visit(
                [](auto& v) -> py::object {
                  return py::cast(*v, py::return_value_policy::copy);
                },
                *tokval_type(x.data()));
          },
          [](TokVal& x, TokVal y) { std::swap(x, y); }, py::doc("The current token value"))
      .def("__repr__",
           [](py::object& x) { return x.attr("value").attr("__repr__")(); })
      .def("__str__",
           [](py::object& x) { return x.attr("value").attr("__str__")(); })
      .def(py::pickle(
          [](const py::object& self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_unique<TokVal>(py::type::of<TokVal>()(t[0]).cast<TokVal>());
          })).doc() = "Custom class holding any value in token form inside ARTS - do not use manually";
  py::implicitly_convertible<WorkspaceVariablesVariant, TokVal>();
  py::implicitly_convertible<py::int_, TokVal>();
  py::implicitly_convertible<py::float_, TokVal>();

  py::class_<MRecord>(m, "MRecord")
      .def_property_readonly("id", &MRecord::Id, ":class:`~pyarts.arts.Index` Method record identity")
      .def_property_readonly("output", &MRecord::Out, ":class:`~pyarts.arts.ArrayOfIndex` Method record output")
      .def_property_readonly("input", &MRecord::In, ":class:`~pyarts.arts.ArrayOfIndex` Method record input")
      .def_property_readonly("value", &MRecord::SetValue, ":class:`~pyarts.arts.TokVal` Method record token values to set")
      .def_property_readonly("tasks", &MRecord::Tasks, ":class:`~pyarts.arts.Agenda` Method record agenda to set")
      .PythonInterfaceBasicRepresentation(MRecord).doc() = "Method record for :class:`~pyarts.arts.Agenda` - do not use manually";

  py::class_<Array<MRecord>>(m, "ArrayOfMRecord")
      .PythonInterfaceArrayDefault(MRecord)
      .PythonInterfaceBasicRepresentation(Array<MRecord>).doc() = ":class:`list` of :class:`~pyarts.arts.MRecord` - do not use manually";

  py::class_<MdRecord>(m, "MdRecord")
      .def(py::init([]() { return std::make_unique<MdRecord>(); }), py::doc("Create empty"))
      .def_property_readonly("name", &MdRecord::Name, ":class:`~pyarts.arts.String` Name of method")
      .def_property_readonly("outs", &MdRecord::Out, ":class:`~pyarts.arts.ArrayOfIndex` Output of method")
      .def_property_readonly("g_out", &MdRecord::GOut, ":class:`~pyarts.arts.ArrayOfString` Generic output of method")
      .def_property_readonly("g_out_types", &MdRecord::GOutType, ":class:`~pyarts.arts.ArrayOfIndex` Types of the generic output")
      .def_property_readonly("ins", &MdRecord::In, ":class:`~pyarts.arts.ArrayOfIndex`  Input to method")
      .def_property_readonly("g_in", &MdRecord::GIn, ":class:`~pyarts.arts.ArrayOfString` Generic input to method")
      .def_property_readonly("g_in_types", &MdRecord::GInType, ":class:`~pyarts.arts.ArrayOfIndex` Types of the generic input")
      .def_property_readonly("g_in_default", &MdRecord::GInDefault, ":class:`~pyarts.arts.ArrayOfString` Generic input parsable defaults to method")
      .def_property_readonly("desc", &MdRecord::Description, ":class:`~pyarts.arts.String` Description of method")
      .PythonInterfaceCopyValue(MdRecord)
      .def("__str__", [](MdRecord& mr) { return var_string('"', mr, '"'); })
      .def("__repr__", [](MdRecord& mr) { return var_string('"', mr, '"'); }).doc() = "Method record for the workspace - do not use manually";

  py::class_<Array<MdRecord>>(m, "ArrayOfMdRecord")
      .PythonInterfaceArrayDefault(MdRecord)
      .PythonInterfaceBasicRepresentation(Array<MdRecord>).doc() = ":class:`list` of :class:`~pyarts.arts.MdRecord` - do not use manually";

  py::class_<AgRecord>(m, "AgRecord")
      .def(py::init([]() { return std::make_unique<AgRecord>(); }), py::doc("Create empty"))
      .PythonInterfaceCopyValue(AgRecord)
      .def_property_readonly("name", &AgRecord::Name, ":class:`~pyarts.arts.String` Name of record")
      .def_property_readonly("description", &AgRecord::Description, ":class:`~pyarts.arts.String` Description of record")
      .def_property_readonly("outs", &AgRecord::Out, ":class:`~pyarts.arts.ArrayOfIndex`  Outputs of record")
      .def_property_readonly("ins", &AgRecord::In, ":class:`~pyarts.arts.ArrayOfIndex`  Inputs to record")
      .def("__str__", [](AgRecord& ar) { return var_string(ar); })
      .def("__repr__", [](AgRecord& ar) { return var_string(ar); }).doc() = "Internal agenda record - do not use manually";

  py::class_<Array<AgRecord>>(m, "ArrayOfAgRecord")
      .PythonInterfaceArrayDefault(AgRecord)
      .PythonInterfaceBasicRepresentation(Array<AgRecord>).doc() = ":class:`list` of :class:`~pyarts.arts.AgRecord` - do not use manually";

  py::class_<Agenda>(m, "Agenda")
      .def(py::init([](Workspace& ws) { return std::make_unique<Agenda>(ws); }), "Create empty")
      .def(py::init([](Workspace&, const Agenda& a) { return a; }),
           py::doc("Copy with extra argument (to mimic :class:`~pyarts.workspace.DelayedAgenda`)"))
      .def(py::init([](Workspace& ws,
                       const std::function<py::object(Workspace&)>& f) {
             return py::cast<Agenda>(f(ws));
           }),
           py::keep_alive<0, 1>(), "Parse :class:`~pyarts.workspace.DelayedAgenda`")
      .def_property_readonly("workspace",
                             py::cpp_function(
                                 [](Agenda& a) {
                                   return a.workspace();
                                 },
                                 py::return_value_policy::reference,
                                 py::keep_alive<0, 1>()),
                             py::doc("Returns an internal workspace - do not use manually"))
      .PythonInterfaceWorkspaceVariableConversion(Agenda)
      .def(py::init([](Workspace& w, const std::filesystem::path& path) {
             return parse_agenda(w,
                                 correct_include_path(path).c_str(),
                                 *w.get<Verbosity>("verbosity"));
           }),
           py::keep_alive<0, 1>(), py::doc("Parse file"))
      .PythonInterfaceFileIO(Agenda)
      .def_property_readonly("main", &Agenda::is_main_agenda, ":class:`bool` Main agenda flag")
      .def_property("name", &Agenda::name, &Agenda::set_name, py::doc(":class:`pyarts.arts.String` Name of agenda"))
      .def_property_readonly("output2push", &Agenda::get_output2push, py::doc(":class:`pyarts.arts.ArrayOfIndex` Additional output"))
      .def_property_readonly("output2dup", &Agenda::get_output2dup, py::doc(":class:`pyarts.arts.ArrayOfIndex` Duplicational output"))
      .def("set_outputs_to_push_and_dup",
           [](Agenda& a) {
             a.set_outputs_to_push_and_dup(
                 *a.workspace()->get<Verbosity>("verbosity"));
           }, py::doc("Sets variable from state"))
      .def(
          "add_workspace_method",
          [](Agenda& a,
             const char* name,
             const py::args& args,
             const py::kwargs& kwargs) {
            auto ws = a.workspace();

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
            Array<AgendaMethodVariable> var_order = sorted_mdrecord(*ws, ptr);

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
                      WorkspaceVariable(*ws, var.group, args[i], var.input).pos;
                } else if (const auto* key = var.name.c_str();
                           kwargs.contains(key)) {
                  var.ws_pos =
                      WorkspaceVariable(*ws, var.group, kwargs[key], var.input)
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
                      error_msg(name, *ws, var_order),
                      "\n"
                      "The generic variable \"",
                      var.name,
                      "\" cannot be converted to a type that matches any available workspace method signatures (note that the signature above is therefore potentially extra weird)")

                  // We have a new variable order but previous variables must be respected
                  auto new_var_order = sorted_mdrecord(*ws, ptr);
                  for (Index j = 0; j < i; j++) {
                    if (var_order[j].ws_pos >= 0 and
                        new_var_order[j].group not_eq
                            WorkspaceVariable(*ws, var_order[j].ws_pos)
                                .group()) {
                      goto increment_ptr;
                    }
                  }

                  // Try with a new type
                  var.group = new_var_order[i].group;
                  goto try_to_set_workspace_position;
                } else {
                  ARTS_USER_ERROR("\nCannot add workspace method:\n",
                                  error_msg(name, *ws, var_order),
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
                                 error_msg(name, *ws, var_order))
            }

            // Adjust pointer to the correct one given the current args and kwargs if we are generic
            if (generic_function) {
              // Ensure we have GIN and GOUT for all of the generic types
              for (Index i = 0; i < var_order.nelem(); i++) {
                auto& var = var_order[i];
                var.group = ws->wsv_data_ptr->operator[](var.ws_pos).Group();
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
                  error_msg(name, *ws, var_order),
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
                        return WorkspaceVariable(*ws, var.ws_pos).group() not_eq
                               var.group;
                      }),
                  "Mismatching signature.  Got:\n",
                  error_msg(name, *ws, var_order),
                  "\n"
                  "But this does not match a signature in Arts")
            }

            // Set the actual values using automatic set methods
            for (auto& var : var_order) {
              WorkspaceVariable val(*ws, var.ws_pos);
              if (std::strncmp(val.name(), "::anon::", 8) == 0) {
                a.push_back(simple_set_method(val));
              }
            }

            // Set the method record
            const auto [in, out] = split_io(var_order);
            const auto m_id = std::distance(global_data::md_data.begin(), ptr);
            if (ptr->SetMethod() and in.nelem()) {
              ARTS_USER_ERROR_IF(
                  ws->is_initialized(in.back()),
                  "Can only set values, use Copy to copy workspace variables")

              // The previous value will be the GIN set value as the last value
              // is also the last in the push_back above
              ARTS_ASSERT(a.Methods().nelem())
              a.push_back(MRecord(m_id,
                                  out,
                                  in,
                                  a.Methods().back().SetValue(),
                                  a.Methods().back().Tasks()));
            } else {
              a.push_back(MRecord(
                  m_id,
                  out,
                  in,
                  {},
                  Agenda{*ws}));
            }

            // Clean up the value
            for (auto& var : var_order) {
              WorkspaceVariable val(*ws, var.ws_pos);
              if (std::strncmp(val.name(), "::anon::", 8) == 0) {
                a.push_back(simple_delete_method(val));
              }
            }
          },
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
          [](Agenda& a, const CallbackFunction& f) mutable {
            auto ws = a.workspace();

            const Index group_index = get_wsv_group_id("CallbackFunction");
            ARTS_USER_ERROR_IF(group_index < 0,
                               "Cannot recognize CallbackFunction")

            const auto counter = ws->nelem();
            const String name = var_string("::callback::", counter);

            Index in = ws->add_wsv(WsvRecord(
                name.c_str(), "Callback created by pybind11 API", group_index));
            ws->push_move(in, std::make_shared<CallbackFunction>(f));

            auto method_ptr =
                global_data::MdMap.find("CallbackFunctionExecute");
            ARTS_USER_ERROR_IF(method_ptr == global_data::MdMap.end(),
                               "Cannot find CallbackFunctionExecute")

            WorkspaceVariable val(*ws, in);
            a.push_back(simple_set_method(val));
            a.push_back(MRecord(
                method_ptr->second,
                {},
                {in},
                {},
                Agenda{val.ws}));
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
          py::arg("f"))
      .def(
          "append_agenda_methods",
          [](Agenda& self, const Agenda& other) {
            ARTS_USER_ERROR_IF(not self.has_same_origin(*other.workspace()), "Agendas on different workspaces")
            const Index n = other.Methods().nelem();  // if other==self
            for (Index i = 0; i < n; i++) self.push_back(other.Methods()[i]);
          },
          py::doc(R"--(
Appends the input agenda's methods to this agenda's method list

Parameters
----------
other : Agenda
    The other agenda

Warning
-------
Both agendas must be defined on the same workspace)--"),
          py::arg("other"))
      .def(
          "execute",
          [](Agenda& a, Workspace& ws) {
            a.set_main_agenda();
            a.execute(ws);
          },
          py::keep_alive<1, 2>(),
          "Executes the agenda as if it was the main agenda")
      .def(
          "check",
          [](Agenda& a, Workspace& w) {
            a.check(w,
                    *static_cast<Verbosity*>(w.get<Verbosity>("verbosity")));
          },
          py::doc("Checks if the agenda works"))
      .def("__repr__",
           [](Agenda& a) { return var_string("Agenda ", a.name()); })
      .def("__str__",
           [](Agenda& a) {
             std::string out =
                 var_string("Arts Agenda ", a.name(), " with methods:\n");
             std::ostringstream os;
             a.print(os, "   ");
             out += os.str();
             return out;
           })
      .def_property("methods", &Agenda::Methods, &Agenda::set_methods, py::doc(":class:`pyarts.arts.ArrayOfMRecord` Methods of the agenda"))
      .def("has_same_origin", &Agenda::has_same_origin, "Checks if another workspace works with this agenda")
      .PythonInterfaceWorkspaceDocumentation(Agenda);

  py::class_<ArrayOfAgenda>(m, "ArrayOfAgenda")
      .def(py::init([]() { return std::make_unique<ArrayOfAgenda>(); }), "Create empty")
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfAgenda)
      .PythonInterfaceCopyValue(ArrayOfAgenda)
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
      }), "Create from :class:`list`")
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
          "check",
          [](ArrayOfAgenda& aa, Workspace& w) {
            for (auto& a : aa)
              a.check(
                  w,
                  *static_cast<Verbosity*>(w.get<Verbosity>("verbosity")));
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
            return std::make_unique<ArrayOfAgenda>(t[0].cast<std::vector<Agenda>>());
          }))
      .PythonInterfaceWorkspaceDocumentationExtra(ArrayOfAgenda, "\n\nThese arrays are partial to contain inter-item logic, please be cautious using them");
  py::implicitly_convertible<std::vector<Agenda>, ArrayOfAgenda>();
}
}  // namespace Python
