#include <auto_md.h>
#include <global_data.h>

#include <algorithm>
#include <fstream>
#include <ostream>
#include <sstream>
#include <utility>
#include <vector>

#include "array.h"
#include "debug.h"
#include "mystring.h"
#include "tokval_io.h"

/** Implements pybind11 bindings generation for Arts workspace
 *
 *  Some important notes:
 *  - The workspace always class definition is always defined as the "ws" variable
 *    - You can define properties and methods onto this class, including read-only and static ones
 *      - The methods are defined using "ws.def(...);"
 *      - The properties are defined using "ws.def_property(...);"
 *      - The first input to such a method is the Workspace as a reference (pointer or true reference)
 *        - Its name cannot be "ws" as well
 *      - The rest of the inputs are what will be seen if calling ws.Method(...) in python
 *  - Use py::doc("Doc") to give documentation to your added method
 *  - Use py::arg("X") or py::arg_v("X", default, "default doc") to give documentation to variables
 *  - Multiple input methods use std::variant
 *    - pybind11 knows how to deal with the input from python, it will find an option
 *    - C++ must deal with selecting the right type
 *    - There are severeal select_* functions in python_interface.h to use, don't hesitate to add more if necessary but try to use them first
 *  - Custom input that can be understood as Arts types use py::object as their type
 *  - Optional inputs (i.e., when there are defaults) make use of std::optional
 *    - We combine std::optional and std::variant when possible
 *  - All workspace variables are of the type WorkspaceVariable
 *  - A workspace variable can be implicitly converted to its base type
 *    - The conversion from python type to C++ prefers the first seen option
 *    - WorkspaceVariable are therefore prefered
 */

struct Group {
  std::string varname_group;
  std::string varname_desc;
  std::size_t artspos;
};

struct Method {
  struct Gin {
    std::vector<std::string> desc;
    std::vector<std::string> group;
    std::vector<std::string> name;
    std::vector<std::string> defs;
    std::vector<bool> hasdefs;
  };

  struct Gout {
    std::vector<std::string> desc;
    std::vector<std::string> group;
    std::vector<std::string> name;
  };

  struct In {
    std::vector<std::string> varname;
    std::vector<std::size_t> varpos;
  };

  struct Out {
    std::vector<std::string> varname;
    std::vector<std::size_t> varpos;
  };

  std::string name;
  std::size_t pos;
  std::string desc;
  std::vector<std::string> authors;
  bool set_method;
  bool agenda_method;
  bool supergeneric;
  bool uses_templates;
  bool pass_workspace;
  bool pass_wsv_names;
  In in;
  Gin gin;
  Out out;
  Gout gout;
  std::vector<std::size_t> inoutvarpos;
};

struct AgendaData {
  std::size_t pos;
  std::string desc;
  std::vector<std::string> ins;
  std::vector<std::string> outs;
};

std::map<std::string, Group> groups() {
  std::map<std::string, std::size_t> group;
  for (auto& x : global_data::WsvGroupMap) group[x.first] = x.second;
  std::map<std::string, std::size_t> name;
  for (auto& x : Workspace::wsv_data) name[x.Name()] = x.Group();
  std::map<std::string, std::string> desc;
  for (auto& x : Workspace::wsv_data) {
    auto& val = desc[x.Name()];
    val = x.Description();
    if (x.has_defaults())
      val += var_string("\nDefault: ", [](auto&& tokval) {
        std::string out = var_string(TokValPrinter{tokval});
        if (out.length() == 0) out += "[]";
        return out;
      }(x.default_value()));
  }
  std::map<std::string, std::size_t> pos;
  for (auto& x : Workspace::WsvMap) pos[x.first] = x.second;

  std::map<std::string, Group> out;
  for (auto& x : name) {
    for (auto& y : group) {
      if (y.second == x.second) {
        out[x.first] = {y.first, desc[x.first], pos[x.first]};
        break;
      }
    }
  }
  return out;
}

std::pair<std::vector<std::string>, std::vector<bool>> fixed_defaults(
    const std::vector<std::string>& vargroups,
    const std::vector<std::string>& vardefaults) {
  std::vector<std::string> defaults(vargroups.size());
  std::vector<bool> hasdefaults(vargroups.size());
  for (size_t i = 0; i < vargroups.size(); i++) {
    if (vardefaults[i] == NODEF)
      hasdefaults[i] = false;
    else
      hasdefaults[i] = true;

    if (vargroups[i] == "String") {
      defaults[i] = std::string("\"") + vardefaults[i] + std::string("\"");
    } else if (vargroups[i] == "Numeric") {
      if ("NaN" == vardefaults[i] or "nan" == vardefaults[i]) {
        defaults[i] = "std::numeric_limits<Numeric>::quiet_NaN()";
      } else if ("Inf" == vardefaults[i] or "inf" == vardefaults[i]) {
        defaults[i] = "std::numeric_limits<Numeric>::infinity()";
      } else if ("-Inf" == vardefaults[i] or "-inf" == vardefaults[i]) {
        defaults[i] = "-std::numeric_limits<Numeric>::infinity()";
      } else {
        defaults[i] = vardefaults[i];
      }
    } else if (vardefaults[i] == "[]") {
      defaults[i] = var_string(vargroups[i], "{}");
    } else {
      defaults[i] = vardefaults[i];
    }

    if (defaults[i] == "") defaults[i] = var_string(vargroups[i], "{}");

    for (auto& x : defaults[i]) {
      if (x == '[')
        x = '{';
      else if (x == ']')
        x = '}';
    }
  }

  return {defaults, hasdefaults};
}

std::vector<Method> methods() {
  std::map<std::string, std::size_t> vargroup;
  for (auto& x : global_data::WsvGroupMap) vargroup[x.first] = x.second;
  std::map<std::string, std::size_t> varpos;
  for (auto& x : Workspace::WsvMap) varpos[x.first] = x.second;
  std::map<std::string, std::size_t> methodpos;
  for (auto& x : global_data::MdMap) methodpos[x.first] = x.second;

  std::vector<std::string> metname;
  for (auto& x : global_data::md_data) metname.push_back(x.Name());
  std::vector<std::string> actual_groups;
  for (auto& x : global_data::md_data)
    actual_groups.push_back(x.ActualGroups());
  std::vector<std::vector<std::size_t>> gin_group;
  for (auto& x : global_data::md_data)
    gin_group.emplace_back(x.GInType().cbegin(), x.GInType().cend());
  std::vector<std::vector<std::string>> gin_names;
  for (auto& x : global_data::md_data)
    gin_names.emplace_back(x.GIn().cbegin(), x.GIn().cend());
  std::vector<std::vector<std::string>> gin_defaults;
  for (auto& x : global_data::md_data)
    gin_defaults.emplace_back(x.GInDefault().cbegin(), x.GInDefault().cend());
  std::vector<std::vector<std::string>> gin_desc;
  for (auto& x : global_data::md_data)
    gin_desc.emplace_back(x.GInDescription().cbegin(),
                          x.GInDescription().cend());
  std::vector<std::vector<std::size_t>> gout_group;
  for (auto& x : global_data::md_data)
    gout_group.emplace_back(x.GOutType().cbegin(), x.GOutType().cend());
  std::vector<std::vector<std::string>> gout_names;
  for (auto& x : global_data::md_data)
    gout_names.emplace_back(x.GOut().cbegin(), x.GOut().cend());
  std::vector<std::vector<std::string>> gout_desc;
  for (auto& x : global_data::md_data)
    gout_desc.emplace_back(x.GOutDescription().cbegin(),
                           x.GOutDescription().cend());
  std::vector<std::vector<std::size_t>> in_wspace;
  for (auto& x : global_data::md_data)
    in_wspace.emplace_back(x.In().cbegin(), x.In().cend());
  std::vector<std::vector<std::size_t>> out_wspace;
  for (auto& x : global_data::md_data)
    out_wspace.emplace_back(x.Out().cbegin(), x.Out().cend());
  std::vector<std::string> desc;
  for (auto& x : global_data::md_data) desc.push_back(x.Description());
  std::vector<std::vector<std::string>> authors;
  for (auto& x : global_data::md_data)
    authors.emplace_back(x.Authors().cbegin(), x.Authors().cend());

  std::vector<bool> set_method;
  for (auto& x : global_data::md_data) set_method.push_back(x.SetMethod());
  std::vector<bool> agenda_method;
  for (auto& x : global_data::md_data)
    agenda_method.push_back(x.AgendaMethod());
  std::vector<bool> supergeneric;
  for (auto& x : global_data::md_data) supergeneric.push_back(x.Supergeneric());
  std::vector<bool> uses_templates;
  for (auto& x : global_data::md_data)
    uses_templates.push_back(x.UsesTemplates());
  std::vector<bool> pass_workspace;
  for (auto& x : global_data::md_data)
    pass_workspace.push_back(x.PassWorkspace());
  std::vector<bool> pass_wsv_names;
  for (auto& x : global_data::md_data)
    pass_wsv_names.push_back(x.PassWsvNames());

  std::vector<std::vector<std::size_t>> inoutvarpos;
  for (auto& x : global_data::md_data)
    inoutvarpos.emplace_back(x.InOut().cbegin(), x.InOut().cend());
  std::vector<std::vector<std::size_t>> invarpos;
  for (auto& x : global_data::md_data)
    invarpos.emplace_back(x.InOnly().cbegin(), x.InOnly().cend());
  std::vector<std::vector<std::size_t>> outvarpos;
  for (auto& x : global_data::md_data)
    outvarpos.emplace_back(x.Out().cbegin(), x.Out().cend());

  std::vector<Method> retval;
  for (std::size_t i = 0; i < desc.size(); i++) {
    Method m;
    m.name = metname[i];
    if (supergeneric[i])
      m.pos = methodpos[metname[i] + String("_sg_") + actual_groups[i]];
    else
      m.pos = methodpos[metname[i]];
    m.desc = desc[i];
    m.authors = authors[i];

    Method::Gin gin;
    gin.desc = gin_desc[i];
    for (auto g : gin_group[i]) {
      bool found = false;
      for (auto& y : vargroup) {
        if (y.second == g) {
          gin.group.push_back(y.first);
          found = true;
          break;
        }
      }
      if (not found) {
        std::cerr << "Cannot find group\n";
        std::terminate();
      }
    }
    gin.name = gin_names[i];
    auto fixgin = fixed_defaults(gin.group, gin_defaults[i]);
    gin.defs = fixgin.first;
    gin.hasdefs = fixgin.second;
    m.gin = gin;

    Method::Gout gout;
    gout.desc = gout_desc[i];
    for (auto g : gout_group[i]) {
      bool found = false;
      for (auto& y : vargroup) {
        if (y.second == g) {
          gout.group.push_back(y.first);
          found = true;
          break;
        }
      }
      if (not found) {
        std::cerr << "Cannot find group\n";
        std::terminate();
      }
    }
    gout.name = gout_names[i];
    m.gout = gout;

    Method::In in;
    for (auto v : in_wspace[i]) {
      bool found = false;
      for (auto& y : varpos) {
        if (v == y.second) {
          in.varname.push_back(y.first);
          found = true;
          break;
        }
      }
      if (not found) {
        std::cerr << "Cannot find variable\n";
        std::terminate();
      }
    }
    in.varpos = invarpos[i];
    m.in = in;

    Method::Out out;
    for (auto v : out_wspace[i]) {
      bool found = false;
      for (auto& y : varpos) {
        if (v == y.second) {
          out.varname.push_back(y.first);
          found = true;
          break;
        }
      }
      if (not found) {
        std::cerr << "Cannot find variable\n";
        std::terminate();
      }
    }
    out.varpos = outvarpos[i];
    m.out = out;

    m.set_method = set_method[i];
    m.agenda_method = agenda_method[i];
    m.supergeneric = supergeneric[i];
    m.uses_templates = uses_templates[i];
    m.pass_workspace = pass_workspace[i];
    m.pass_wsv_names = pass_wsv_names[i];
    m.inoutvarpos = inoutvarpos[i];
    retval.push_back(m);
  }

  return retval;
}

std::map<std::string, AgendaData> agendas() {
  auto g = groups();

  std::map<std::string, AgendaData> out;
  for (auto& x : global_data::agenda_data) {
    out[x.Name()].pos = global_data::AgendaMap.at(x.Name());
    out[x.Name()].desc = x.Description();

    for (std::size_t i : x.In()) {
      bool found = false;
      for (auto& y : g) {
        if (y.second.artspos == i) {
          out[x.Name()].ins.push_back(y.first);
          found = true;
          break;
        }
      }
      if (not found) {
        std::cerr << "Cannot find the variable\n";
        std::terminate();
      }
    }

    for (std::size_t i : x.Out()) {
      bool found = false;
      for (auto& y : g) {
        if (y.second.artspos == i) {
          out[x.Name()].outs.push_back(y.first);
          found = true;
          break;
        }
      }
      if (not found) {
        std::cerr << "Cannot find the variable\n";
        std::terminate();
      }
    }
  }
  return out;
}

struct NameMaps {
  std::map<std::string, AgendaData> agendaname_agenda;
  std::vector<Method> methodname_method;
  std::map<std::string, Group> varname_group;
  std::map<std::string, std::size_t> group;

  NameMaps() {
    for (auto& x : global_data::WsvGroupMap) group[x.first] = x.second;
    varname_group = groups();
    methodname_method = methods();
    agendaname_agenda = agendas();
  }
};

void includes(std::ofstream& os) {
  os << "#include <python_interface.h>" << '\n';
}

void workspace_variables(size_t n, const NameMaps& arts) {
  const Index verbpos = arts.varname_group.find("verbosity")->second.artspos;

  std::vector<std::ofstream> oss(n);
  for (size_t i = 0; i < n; i++) {
    oss[i] =
        std::ofstream(var_string("py_auto_workspace_split_vars_", i, ".cc"));
    oss[i]
        << "#include <python_interface.h>\n#include <pybind11/functional.h>\n\nnamespace Python {\nvoid py_auto_workspace_wsv_"
        << i << "(py::class_<Workspace>& ws [[maybe_unused]]) {\n";
  }

  auto osptr = oss.begin();
  for (auto& [name, data] : arts.varname_group) {
    auto& os = *osptr;

    os << "ws.def_property(\"" << name << "\",\n";
    os << R"--(  py::cpp_function([](Workspace& w) -> WorkspaceVariable {
    return WorkspaceVariable{w, )--"
       << data.artspos
       << "};\n  }, py::return_value_policy::reference_internal),\n";

    os << "  [](Workspace& w, " << data.varname_group << "& val) {\n";

    if (data.varname_group == "Agenda")
      os << "    val.set_name(\"" << name
         << "\");\n    val.check(w, *static_cast<Verbosity*>(w[" << verbpos
         << "].get()));\n";
    if (data.varname_group == "ArrayOfAgenda")
      os << "    for(auto& a : val) a.set_name(\"" << name
         << "\");\n    for(auto& a : val) a.check(w, *static_cast<Verbosity*>(w["
         << verbpos << "].get()));\n";
    os << "    " << data.varname_group << "& v = WorkspaceVariable {w, "
       << data.artspos << "};\n    v = val;\n  }, py::doc(R\"-x-("
       << data.varname_desc << ")-x-\")\n); \n\n";

    osptr++;
    if (osptr == oss.end()) osptr = oss.begin();
  }

  for (auto& os : oss) os << "}\n}\n\n";
}

void print_method_desc(std::ofstream& os,
                       const Method& method,
                       const std::map<std::string, Group>& groups,
                       bool pass_verbosity) {
  os << ",\npy::doc(\nR\"-METHODS_DESC-(\n" << method.desc;
  os << "\nAuthor(s): ";
  for (auto& author : method.authors) {
    if (&author == &method.authors.back() and method.authors.size() > 1) {
      if (method.authors.size() > 2) os << ',';
      os << " and ";
    } else if (&author not_eq &method.authors.front())
      os << ", ";
    os << author;
  }
  os << '\n' << '\n' << "Parameters\n----------\n";

  for (const auto& i : method.out.varname) {
    os << i << " : " << "pyarts.arts." << groups.at(i).varname_group << ", optional\n";
    os << "    As WSV (";
    if (std::none_of(method.in.varname.cbegin(),
                     method.in.varname.cend(),
                     [out = i](const auto& in) { return in == out; })) {
      os << "IN";
    }
    os << "OUT)\n";
  }
  for (size_t i = 0; i < method.gout.name.size(); i++) {
    os << method.gout.name[i] << " : " << "pyarts.arts." << method.gout.group[i] << "\n    "
       << method.gout.desc[i] << " (OUT)\n";
  }
  for (const auto& i : method.in.varname) {
    if (std::none_of(method.out.varname.cbegin(),
                     method.out.varname.cend(),
                     [in = i](const auto& out) { return in == out; })) {
      os << i << " : " << "pyarts.arts." << groups.at(i).varname_group
         << ", optional\n    As WSV (IN)\n";
    }
  }
  for (size_t i = 0; i < method.gin.name.size(); i++) {
    os << method.gin.name[i] << " : " << "pyarts.arts." << method.gin.group[i];
    if (method.gin.hasdefs[i]) {
      os << ", optional";
    }
    os << "\n    " << method.gin.desc[i] << " (IN";
    if (method.gin.hasdefs[i]) {
      os << "; default: " << method.gin.defs[i];
    }
    os << ')' << '\n';
  }

  if (pass_verbosity)
    os << "verbosity : pyarts.arts.Verbosity, optional\n    As WSV (IN)\n";

  os << "\n)-METHODS_DESC-\")";
}

void print_method_args(std::ofstream& os,
                       const Method& method,
                       bool pass_verbosity) {
  for (const auto& i : method.out.varname) {
    os << ',' << '\n'
       << "py::arg_v(\"" << i << "\", std::nullopt, \"ws." << i
       << "\").noconvert()";
  }
  for (const auto& i : method.gout.name) {
    os << ',' << '\n' << "py::arg(\"" << i << "\").noconvert().none(false)";
  }
  for (const auto& i : method.in.varname) {
    if (std::none_of(method.out.varname.cbegin(),
                     method.out.varname.cend(),
                     [in = i](const auto& out) { return in == out; })) {
      os << ',' << '\n'
         << "py::arg_v(\"" << i << "\", std::nullopt, \"ws." << i << "\")";
    }
  }
  for (size_t i = 0; i < method.gin.name.size(); i++) {
    if (method.gin.hasdefs[i]) {
      if (method.gin.group[i] == "String")
        os << ',' << '\n'
           << "py::arg_v(\"" << method.gin.name[i] << R"(", std::nullopt, "\)"
           << method.gin.defs[i] << R"("\"")" << ')';
      else
        os << ',' << '\n'
           << "py::arg_v(\"" << method.gin.name[i] << "\", std::nullopt, \""
           << method.gin.defs[i] << "\")";
    } else {
      os << ',' << '\n'
         << "py::arg(\"" << method.gin.name[i] << "\").none(false)";
    }
  }
  if (pass_verbosity)
    os << ',' << '\n'
       << R"--(py::arg_v("verbosity", std::nullopt, "ws.verbosity"))--";
}

void workspace_method_create(size_t n, const NameMaps& arts) {
  std::vector<std::ofstream> oss(n);
  for (size_t i = 0; i < n; i++) {
    oss[i] =
        std::ofstream(var_string("py_auto_workspace_split_create_", i, ".cc"));
    oss[i]
        << "#include \"py_auto_interface.h\"\n#include <pybind11/stl.h>\n\nnamespace Python {\nvoid py_auto_workspace_wsc_"
        << i << "(py::class_<Workspace>& ws [[maybe_unused]]) {\n";
  }
  auto osptr = oss.begin();

  for (auto& [group, group_index] : arts.group) {
    auto& os = *osptr;

    if (group == "Any") continue;
    os << "ws.def(\"" << group
       << "Create\",[](Workspace& w, const char * name, std::optional<const char *> desc, std::optional<"
       << group;
    if (group == "Index" or group == "Numeric") os << '_';
    os << R"--(> value) {
  Index pos=-1;
  if (auto varptr = w.WsvMap.find(name); varptr == w.WsvMap.end()) {
    ARTS_USER_ERROR_IF(std::string_view sv{name}; sv.size() == 0 or sv.front() == ':' or sv.back() == ':',
      "Cannot create variables without a name, or one that starts or ends with \":\"\nVariable name: \"", name, "\"")
    pos = w.add_wsv_inplace(WsvRecord(name, desc.has_value() ? *desc : "User-generated workspace variable", )--"
       << group_index << R"--());
  } else {
    pos = varptr->second;
    ARTS_USER_ERROR_IF(w.wsv_data[pos].Group() not_eq )--"
       << group_index << R"--(, "Already exist of different group: ", name)
  }

  if (value.has_value()) w.push_move(pos, std::make_shared<)--"
       << group << R"--(>(*value));
}, py::doc(R"-x-(
Create new )--"
       << group << R"--( on the workspace

It is recommended that this is only called once per Workspace instance.

If there is no default value, the variable remains uninitialized

If the variable already exist, an error is thrown only if the group
is not )--"
       << group << R"--(.  If the variable does exist, a new value
is pushed onto its stack iff value is given.  This efficiently puts
any current value out of reach for the current Agenda level.  Note
that this should be OK inside Agendas as the stack should be freed
at the end of the Agenda call, however if this is called directly
on the workspace outside an Agenda, the extra memory allocated can
only be freed using "del name" as an appropriate pair

Parameters:
-----------
name (str): Name of the variable, can be accessed later with getattr() on the workspace
desc (str): Description of the variable [Optional]
value ()--"
       << group << R"--(): Initialized value [Optional]
)-x-"),
  py::arg("name").none(false),
  py::arg("desc")=std::nullopt,
  py::arg("value")=std::nullopt);
  
)--";

    osptr++;
    if (osptr == oss.end()) osptr = oss.begin();
  }

  for (auto& os : oss) os << "}\n}\n\n";
}

void workspace_method_nongenerics(size_t n, const NameMaps& arts) {
  std::vector<std::ofstream> oss(n);
  for (size_t i = 0; i < n; i++) {
    oss[i] =
        std::ofstream(var_string("py_auto_workspace_split_methods_", i, ".cc"));
    oss[i]
        << "#include <python_interface.h>\n\nnamespace Python {\nvoid py_auto_workspace_wsm_"
        << i << "(py::class_<Workspace>& ws [[maybe_unused]]) {\n";
  }

  auto osptr = oss.begin();

  for (auto& method : arts.methodname_method) {
    auto& os = *osptr;

    if (method.agenda_method) continue;  // FIXME: Should be fixed???

    // Skip create methods
    if (std::any_of(arts.group.cbegin(),
                    arts.group.cend(),
                    [metname = method.name](auto& y) {
                      return (y.first + String("Create")) == metname;
                    }))
      continue;

    const bool pass_verbosity = std::none_of(
        method.out.varname.begin(), method.out.varname.end(), [](auto& var) {
          return var == "verbosity";
        });

    const bool pass_workspace =
        method.pass_workspace or method.agenda_method or
        std::any_of(
            method.gin.group.cbegin(),
            method.gin.group.cend(),
            [](auto& g) { return g == "Agenda" or g == "ArrayOfAgenda"; }) or
        std::any_of(
            method.in.varname.cbegin(), method.in.varname.cend(), [&](auto& g) {
              return arts.varname_group.at(g).varname_group == "Agenda" or
                     arts.varname_group.at(g).varname_group == "ArrayOfAgenda";
            });

    if ([&]() {
          long i = 0;
          for (auto& a : arts.methodname_method) i += (a.name == method.name);
          return i > 1;
        }())
      continue;  // Skip 'generic' functions

    // Arguments from python side
    os << "ws.def(\"" << method.name << "\",[](Workspace& w_ [[maybe_unused]]";
    for (const auto& i : method.out.varname) {
      auto& group = arts.varname_group.at(i).varname_group;
      os << ",\nstd::optional<std::variant<";
      if (std::any_of(method.in.varname.cbegin(),
                      method.in.varname.cend(),
                      [out = i](const auto& in) { return in == out; }))
        os << "const ";
      os << "WorkspaceVariable *, " << group;
      if (group == "Index" or group == "Numeric") os << "_";
      os << " *>> " << i;
    }
    for (std::size_t i = 0; i < method.gout.name.size(); i++) {
      auto& group = method.gout.group[i];
      os << ",\nstd::variant<WorkspaceVariable *, " << group;
      if (group == "Index" or group == "Numeric") os << "_";
      os << " *> " << method.gout.name[i];
    }
    for (const auto& i : method.in.varname) {
      if (std::none_of(method.out.varname.cbegin(),
                       method.out.varname.cend(),
                       [in = i](const auto& out) { return in == out; })) {
        auto& group = arts.varname_group.at(i).varname_group;
        os << ",\nstd::optional<std::variant<const WorkspaceVariable *, "
           << group;
        if (group == "Index" or group == "Numeric") os << "_";
        os << " *>> " << i;
      }
    }
    for (std::size_t i = 0; i < method.gin.name.size(); i++) {
      auto& group = method.gin.group[i];
      os << ",\n";
      if (method.gin.hasdefs[i]) os << "std::optional<";
      os << "std::variant<const WorkspaceVariable *, " << group;
      if (group == "Index" or group == "Numeric") os << "_";
      os << " *>";
      if (method.gin.hasdefs[i]) os << '>';
      os << ' ' << method.gin.name[i];
    }
    if (pass_verbosity)
      os << ",\nstd::optional<std::variant<const WorkspaceVariable *, Verbosity *>> verbosity";
    os << ") {\n";
    bool has_any = false;

    // Create defaults
    has_any = false;
    for (std::size_t i = 0; i < method.gin.name.size(); i++) {
      if (method.gin.hasdefs[i]) {
        os << "  const auto " << method.gin.name[i]
           << "_ = " << method.gin.group[i] << '(' << method.gin.defs[i] << ')'
           << ";\n";
        has_any = true;
      }
    }
    if (has_any) os << '\n';

    // Outputs
    Index counter = 0;
    for (const auto& i : method.out.varname) {
      auto& group = arts.varname_group.at(i).varname_group;
      os << "  auto& arg" << counter++ << "_ = select_";
      if (std::any_of(method.in.varname.cbegin(),
                      method.in.varname.cend(),
                      [out = i](const auto& in) { return in == out; }))
        os << "in";
      os << "out<";
      os << group;
      os << ">(WorkspaceVariable{w_, " << arts.varname_group.at(i).artspos
         << "}, " << i << ");\n";
    }

    // Generic Outputs
    for (std::size_t i = 0; i < method.gout.name.size(); i++) {
      os << "  auto& arg" << counter++ << "_ = select_gout<";
      os << method.gout.group[i];
      os << ">(" << method.gout.name[i] << ");\n";
    }

    // Generic Output Names
    if (method.pass_wsv_names) {
      for (auto& name : method.gout.name) {
        os << "  const String arg" << counter++
           << "_{std::holds_alternative<WorkspaceVariable*>(" << name
           << ") ? std::get<WorkspaceVariable*>(" << name << ") -> name() : \""
           << name << "\"};\n";
      }
    }

    // Inputs (that are not outputs)
    for (auto& i : method.in.varname) {
      if (std::none_of(method.out.varname.cbegin(),
                       method.out.varname.cend(),
                       [in = i](const auto& out) { return in == out; })) {
        auto& group = arts.varname_group.at(i).varname_group;
        os << "  const auto& arg" << counter++ << "_ = select_in<";
        os << group;
        os << ">(WorkspaceVariable{w_, " << arts.varname_group.at(i).artspos
           << "}, " << i << ");\n";
      }
    }

    // Generic Inputs
    for (std::size_t i = 0; i < method.gin.name.size(); i++) {
      os << "  const auto& arg" << counter++ << "_ = select_gin<"
         << method.gin.group[i] << ">(";
      if (method.gin.hasdefs[i]) os << method.gin.name[i] << "_, ";
      os << method.gin.name[i] << ");\n";
    }

    // Generic Input Names
    if (method.pass_wsv_names) {
      for (std::size_t i = 0; i < method.gin.name.size(); i++) {
        auto& name = method.gin.name[i];
        os << "  const String arg" << counter++ << "_{";
        if (method.gin.hasdefs[i]) os << name << ".has_value() ? (";
        os << "std::holds_alternative<const WorkspaceVariable*>(";
        if (method.gin.hasdefs[i]) os << '*';
        os << name << ") ? std::get<const WorkspaceVariable *>(";
        if (method.gin.hasdefs[i]) os << '*';
        os << name << ") -> name() : \"" << name << "\"";
        if (method.gin.hasdefs[i]) os << ") : \"" << name << "\"";
        os << "};\n";
      }
    }

    // Verbosity?
    if (pass_verbosity) {
      os << "  const auto& arg" << counter++
         << "_ = select_in<Verbosity>(WorkspaceVariable{w_, "
         << arts.varname_group.at("verbosity").artspos << "}, verbosity);\n";
    }
    os << '\n';

    // Arguments from Arts side
    has_any = false;
    os << "  " << method.name << '(';
    if (pass_workspace) {
      os << "w_";
      has_any = true;
    }
    for (Index current_count = 0; current_count < counter; current_count++) {
      if (has_any) os << ", ";
      has_any = true;
      os << "arg" << current_count << "_";
    }
    os << ");\n}";

    // Name the paramters and show their defaults
    print_method_args(os, method, pass_verbosity);

    // Put description at the end
    print_method_desc(os, method, arts.varname_group, pass_verbosity);

    os << ')' << ';' << '\n' << '\n';

    osptr++;
    if (osptr == oss.end()) osptr = oss.begin();
  }

  for (auto& os : oss) os << "}\n}\n\n";
}

struct ArgumentHelper {
  bool in;
  bool out;
  bool gen;
  bool def;
  bool wsv;
  Index wsvpos;
  String name;
  String def_str;
  ArrayOfString types;
};

struct VariadicInfo {
  Index counter;
  String name, type, first, second;
};

void workspace_method_generics(size_t n, const NameMaps& arts) {
  std::vector<std::ofstream> oss(n);
  for (size_t i = 0; i < n; i++) {
    oss[i] =
        std::ofstream(var_string("py_auto_workspace_split_generic_", i, ".cc"));
    oss[i]
        << "#include <python_interface.h>\n\nnamespace Python {\nvoid py_auto_workspace_wsg_"
        << i << "(py::class_<Workspace>& ws [[maybe_unused]]) {\n";
  }

  auto osptr = oss.begin();

  ArrayOfString ignore_methods{"Delete"};

  ArrayOfString extra_workspace_for_agenda{"Append", "Copy"};

  for (std::size_t imethod = 0; imethod < arts.methodname_method.size();
       imethod++) {
    auto& method = arts.methodname_method[imethod];
    auto& os = *osptr;

    if (method.agenda_method) {
      continue;  // FIXME: Should be fixed???
    }

    // Skip create methods
    if (std::any_of(arts.group.cbegin(),
                    arts.group.cend(),
                    [metname = method.name](auto& y) {
                      return (y.first + String("Create")) == metname;
                    }))
      continue;

    if (ignore_methods.end() not_eq
        std::find(ignore_methods.begin(), ignore_methods.end(), method.name))
      continue;

    const bool pass_verbosity = std::none_of(
        method.out.varname.begin(), method.out.varname.end(), [](auto& var) {
          return var == "verbosity";
        });

    const bool pass_workspace =
        method.pass_workspace or method.agenda_method or
        std::any_of(
            method.gin.group.cbegin(),
            method.gin.group.cend(),
            [](auto& g) { return g == "Agenda" or g == "ArrayOfAgenda"; }) or
        std::any_of(
            method.in.varname.cbegin(), method.in.varname.cend(), [&](auto& g) {
              return arts.varname_group.at(g).varname_group == "Agenda" or
                     arts.varname_group.at(g).varname_group == "ArrayOfAgenda";
            });

    // Count the number of generic functions
    const std::size_t n_func = [&]() {
      std::size_t i = 0;
      for (auto& a : arts.methodname_method) i += (a.name == method.name);
      return i;
    }();
    if (n_func < 2) continue;  // Skip 'nongeneric' functions
    const std::size_t n_method_end = imethod + n_func;

    // Fix the argument order already so generic variables can be known
    Array<ArgumentHelper> arg_help;
    bool first = true;
    for (; imethod < n_method_end; imethod++) {
      auto& meth = arts.methodname_method[imethod];
      Index pos = 0;

      for (std::size_t i = 0; i < meth.out.varname.size(); i++) {
        if (first) {
          auto& arg = arg_help.emplace_back();
          arg.in = std::any_of(meth.in.varname.cbegin(),
                               meth.in.varname.cend(),
                               [out = meth.out.varname[i]](const auto& in) {
                                 return in == out;
                               });
          arg.out = true;
          arg.gen = false;
          arg.def = false;
          arg.wsv = true;
          arg.name = meth.out.varname[i];
          arg.wsvpos = arts.varname_group.at(arg.name).artspos;
          arg.def_str = "";
          arg.types =
              ArrayOfString{arts.varname_group.at(arg.name).varname_group};
        }

        pos++;
      }

      for (std::size_t i = 0; i < meth.gout.name.size(); i++) {
        if (first) {
          auto& arg = arg_help.emplace_back();
          arg.in = false;
          arg.out = true;
          arg.gen = true;
          arg.def = false;
          arg.wsv = false;
          arg.wsvpos = -1;
          arg.name = meth.gout.name[i];
          arg.def_str = "";
          arg.types = ArrayOfString(0);
        }
        arg_help[pos].types.push_back(meth.gout.group[i]);

        pos++;
      }

      for (std::size_t i = 0; i < meth.in.varname.size(); i++) {
        if (std::none_of(meth.out.varname.cbegin(),
                         meth.out.varname.cend(),
                         [in = meth.in.varname[i]](const auto& out) {
                           return in == out;
                         })) {
          if (first) {
            auto& arg = arg_help.emplace_back();
            arg.in = true;
            arg.out = false;
            arg.gen = false;
            arg.def = false;
            arg.wsv = true;
            arg.name = meth.in.varname[i];
            arg.wsvpos = arts.varname_group.at(arg.name).artspos;
            arg.def_str = "";
            arg.types =
                ArrayOfString{arts.varname_group.at(arg.name).varname_group};
          }

          pos++;
        }
      }

      for (std::size_t i = 0; i < meth.gin.name.size(); i++) {
        if (first) {
          auto& arg = arg_help.emplace_back();
          arg.in = true;
          arg.out = false;
          arg.gen = true;
          arg.def = meth.gin.hasdefs[i];
          arg.wsv = false;
          arg.wsvpos = -1;
          arg.name = meth.gin.name[i];
          arg.def_str = meth.gin.defs[i];
          arg.types = ArrayOfString(0);
        }
        arg_help[pos].types.push_back(meth.gin.group[i]);

        pos++;
      }
      first = false;
    }
    imethod--;

    // Minimizie generic inputs
    std::for_each(arg_help.begin(), arg_help.end(), [](auto& arg) {
      if (std::all_of(arg.types.begin(),
                      arg.types.end(),
                      [in = arg.types.front()](auto& current) {
                        return current == in;
                      }))
        arg.types = ArrayOfString{arg.types.front()};
    });

    bool allow_py_object_input = false;
    {
      bool has_one_possible_combination = true;
      Index generic_var_out_count = 0;
      Index generic_var_in_count = 0;
      for (auto& arg : arg_help) {
        ArrayOfString copies = arg.types;
        std::sort(copies.begin(), copies.end());
        if (std::adjacent_find(copies.begin(), copies.end()) not_eq
            copies.end())
          has_one_possible_combination = false;
        generic_var_in_count += arg.types.nelem() > 1 and arg.in;
        generic_var_out_count += arg.types.nelem() > 1 and arg.out;
      }

      if (has_one_possible_combination and generic_var_in_count == 1 and
          generic_var_out_count > 0)
        allow_py_object_input = true;
    }

    bool has_any = false;

    // Arguments from python side
    os << "ws.def(\"" << method.name << "\",[](Workspace& w_ [[maybe_unused]]";
    for (auto& arg : arg_help) {
      os << ",\n";
      if (arg.wsv or arg.def) os << "std::optional<";
      os << "std::variant<";
      if (arg.in) os << "const ";
      os << "WorkspaceVariable *, ";
      if (arg.types.size() == 1) {
        os << arg.types.front();
        if (arg.types.front() == "Index" or arg.types.front() == "Numeric")
          os << '_';
        os << " *";
      } else {
        ArrayOfString copy_types = arg.types;
        std::sort(copy_types.begin(), copy_types.end());
        copy_types.erase(std::unique(copy_types.begin(), copy_types.end()),
                         copy_types.end());

        has_any = false;
        for (auto& t : copy_types) {
          if (has_any) os << ", ";
          has_any = true;
          os << t;
          if (t == "Index" or t == "Numeric") os << "_";
          os << " *";
        }

        if (arg.in and not arg.out and allow_py_object_input)
          os << ", py::object *";
      }
      os << '>';

      if (not arg.gen or arg.def) os << '>';
      os << ' ' << arg.name;
    }
    if (pass_verbosity)
      os << ",\nstd::optional<std::variant<const WorkspaceVariable *, Verbosity *>> verbosity";
    os << ") {\n";

    // Create defaults
    has_any = false;
    for (std::size_t i = 0; i < method.gin.name.size(); i++) {
      if (method.gin.hasdefs[i]) {
        os << "  const auto " << method.gin.name[i] << '_' << '='
           << method.gin.group[i] << '{';
        if (method.gin.defs[i] not_eq "{}") os << method.gin.defs[i];
        os << "};\n";
        has_any = true;
      }
    }
    if (has_any) os << '\n';

    // Variables
    Index counter = 0;
    for (auto& arg : arg_help) {
      if (method.pass_wsv_names and arg.gen) {
        os << "  const String arg" << counter << "_name_{";
        if (arg.def) os << arg.name << ".has_value() ? (";
        os << "std::holds_alternative<";
        if (arg.in) os << "const ";
        os << "WorkspaceVariable*>(";
        if (arg.def) os << '*';
        os << arg.name << ") ? std::get<";
        if (arg.in) os << "const ";
        os << "WorkspaceVariable*>(";
        if (arg.def) os << '*';
        os << arg.name << ") -> name() : \"" << arg.name << "\"";
        if (arg.def) os << ") : \"" << arg.name << "\"";
        os << "};\n";
      }

      os << "  ";
      if (arg.types.size() == 1) {
        if (arg.in and not arg.out) os << "const ";
        os << "auto& ";
      } else {
        os << "auto wvv_";
      }
      os << "arg" << counter++ << "_ = select_";
      if (arg.types.size() == 1) {
        if (arg.gen) os << 'g';
        if (arg.in) os << "in";
        if (arg.out) os << "out";
        os << '<' << arg.types.front() << ">(";
        if (arg.wsv) os << "WorkspaceVariable{w_, " << arg.wsvpos << "}, ";
        if (arg.def) os << arg.name << "_, ";
        os << arg.name << ");";
      } else {
        os << "wvv(" << arg.name << ");";
      }
      os << '\n';
    }

    if (pass_verbosity) {
      os << "  const auto& arg" << counter++
         << "_ = select_in<Verbosity>(WorkspaceVariable{w_, "
         << arts.varname_group.at("verbosity").artspos << "}, verbosity);\n";
    }
    os << '\n';

    // Find the right function
    os << "  ";
    for (std::size_t i = 0; i < n_func; i++) {
      os << "if (";
      counter = 0;
      has_any = false;

      for (auto& arg : arg_help) {
        if (arg.types.size() > 1) {
          if (arg.in and not arg.out and allow_py_object_input) continue;

          if (has_any) os << " and ";
          has_any = true;
          os << "std::holds_alternative<" << arg.types[i];
          if (arg.types[i] == "Numeric" or arg.types[i] == "Index") os << '_';
          os << " *>(wvv_arg" << counter << "_)";
        }
        counter++;
      }
      os << ") {\n";

      // Interpret the right variables
      counter = 0;
      std::vector<VariadicInfo> input_var_args;
      for (auto& arg : arg_help) {
        if (arg.types.size() > 1) {
          if (arg.in and not arg.out and allow_py_object_input) {
            VariadicInfo x;
            x.counter = counter;
            x.name = var_string("arg", counter, '_');
            x.type = arg.types[i];
            bool num = (arg.types[i] == "Numeric" or arg.types[i] == "Index");
            x.first = var_string("* _tmp",
                                 counter,
                                 ".cast<",
                                 arg.types[i],
                                 num ? '_' : ' ',
                                 "*>()");
            x.second = var_string("*std::get<",
                                  arg.types[i],
                                  (num ? '_' : ' '),
                                  "*>(wvv_arg",
                                  counter,
                                  "_)");
            input_var_args.push_back(x);
            continue;
          }

          os << "    ";
          if (not arg.out) os << "const ";
          os << arg.types[i];
          os << "& arg" << counter << "_ = ";
          if (arg.in and not arg.out and allow_py_object_input) {
            os << "std::holds_alternative<py::object *>(wvv_arg" << counter
               << "_) ? *std::get<py::object *>(wvv_arg" << counter
               << "_)->cast<" << arg.types[i];
            if (arg.types[i] == "Numeric" or arg.types[i] == "Index") os << '_';
            os << " *>() : ";
          }
          os << "**std::get_if<" << arg.types[i];
          if (arg.types[i] == "Numeric" or arg.types[i] == "Index") os << '_';
          os << " *>(&wvv_arg" << counter << "_);\n";
        }
        counter++;
      }

      // Arguments from Arts side
      has_any = false;
      std::ostringstream method_os;
      method_os << method.name << '(';
      if (pass_workspace or
          (std::any_of(extra_workspace_for_agenda.begin(),
                       extra_workspace_for_agenda.end(),
                       [name = method.name](auto& special_name) {
                         return special_name == name;
                       }) and
           std::any_of(arg_help.begin(), arg_help.end(), [i](auto& arg) {
             return arg.types.size() > 1 and (arg.types[i] == "Agenda" or
                                              arg.types[i] == "ArrayOfAgenda");
           }))) {
        method_os << "w_";
        has_any = true;
      }

      counter = 0;
      Index gins = 0, gouts = 0;
      for (auto& arg : arg_help) {
        gins += arg.gen and arg.in;
        gouts += arg.gen and arg.out;

        while (method.pass_wsv_names and not arg.out and gouts not_eq 0) {
          method_os << ", arg" << counter - (gouts--) << "_name_";
        }

        if (has_any) method_os << ", ";
        has_any = true;
        method_os << "arg" << counter++ << "_";
      }

      while (method.pass_wsv_names and gins not_eq 0) {
        method_os << ", arg" << counter - (gins--) << "_name_";
      }

      if (pass_verbosity) method_os << ", arg" << counter << "_";
      method_os << ");";

      const String method_call = method_os.str();

      if (allow_py_object_input and input_var_args.size()) {
        Index s = 4;
        auto spaces = [](Index NUM) {
          std::ostringstream space;
          while (NUM--) space << ' ';
          return space.str();
        };

        // keywords to find a logical path
        const String path1 = "__#1__";
        const String path2 = "__#2__";
        String path1_exe;
        String path2_exe;

        // Get the last argument to work from
        VariadicInfo last = input_var_args.back();
        input_var_args.pop_back();
        ARTS_ASSERT(
            input_var_args.size() == 0,
            "You need to implement a better parser. Sorry!  We had no examples of multiple "
            "GIN in a function that gives GOUT before you added the ",
            method.name,
            " method, so this is not supported (yet) as it could not be tested")

        String full_exe =
            var_string(spaces(s),
                       "if (std::holds_alternative<py::object *>(wvv_arg",
                       last.counter,
                       "_)) {\n",
                       spaces(s + 2),
                       path1,
                       '\n',
                       spaces(s),
                       "} else {\n",
                       spaces(s + 2),
                       path2,
                       '\n',
                       spaces(s),
                       '}',
                       '\n');

        // Path1 may have to be formed if it contains an agenda
        if (last.type == "Agenda") {
          path1_exe = var_string(
              "auto _tmp",
              last.counter,
              " = py::type::of<Agenda>()(w_, * std::get<py::object *>(wvv_arg",
              last.counter,
              "_));\n",
              spaces(s + 2),
              "const Agenda& ",
              last.name,
              " = * _tmp",
              last.counter,
              ".cast<Agenda *>();\n",
              spaces(s + 2),
              method_call);
        } else {
          bool num = (last.type == "Numeric" or last.type == "Index");
          path1_exe = var_string("auto _tmp",
                                 last.counter,
                                 " = py::type::of<",
                                 last.type,
                                 num ? '_' : ' ',
                                 ">()(* std::get<py::object *>(wvv_arg",
                                 last.counter,
                                 "_));\n",
                                 spaces(s + 2),
                                 "const ",
                                 last.type,
                                 '&',
                                 ' ',
                                 last.name,
                                 " = ",
                                 last.first,
                                 ';',
                                 '\n',
                                 spaces(s + 2),
                                 method_call);
        }

        // Make sure we can get references out of these values
        path2_exe = var_string("const ",
                               last.type,
                               "& ",
                               last.name,
                               " = ",
                               last.second,
                               ';',
                               '\n',
                               spaces(s + 2),
                               method_call);

        full_exe.replace(full_exe.find(path1), path1.length(), path1_exe);
        full_exe.replace(full_exe.find(path2), path2.length(), path2_exe);
        //std::cerr << full_exe << '\n';

        os << full_exe;
      } else {
        os << "    " << method_call << '\n';
      }
      os << "  } else ";
    }

    os << "{\n"
          "    WorkspaceVariablesVariant2String printer{};\n"
          "    ARTS_USER_ERROR(\"Call to invalid signature:\"\n"
          "                    \""
       << method.name << '(';
    has_any = false;
    counter = 0;
    for (auto& arg : arg_help) {
      if (has_any) os << ", ";
      has_any = true;
      os << arg.name << " : ";
      if (arg.types.size() > 1)
        os << "\",\n"
              "                    "
              "std::visit(printer, wvv_arg"
           << counter
           << "_),\n"
              "                    \"";
      else
        os << arg.types.front();
      counter++;
    }
    if (pass_verbosity) {
      if (has_any) os << ", ";
      os << "verbosity : Verbosity";
    }
    os << ")\"\n"
          "                    \"\\n\\n\"\n"
          "                    \"See method description for valid signatures\")\n  }\n}";

    // Name the paramters and show their defaults
    print_method_args(os, method, pass_verbosity);

    // Put description at the end
    print_method_desc(os, method, arts.varname_group, pass_verbosity);

    os << ')' << ';' << '\n' << '\n';

    osptr++;
    if (osptr == oss.end()) osptr = oss.begin();
  }

  // Delete
  {
    auto& os = *osptr;
    const Index verbpos = arts.varname_group.find("verbosity")->second.artspos;
    auto method_ptr = std::find_if(arts.methodname_method.begin(),
                                   arts.methodname_method.end(),
                                   [](auto& m) { return m.name == "Delete"; });
    const bool pass_verbosity =
        std::none_of(method_ptr->out.varname.begin(),
                     method_ptr->out.varname.end(),
                     [](auto& var) { return var == "verbosity"; });
    ARTS_USER_ERROR_IF(method_ptr == arts.methodname_method.end(),
                       "The Delete method no longer exist")

    os << R"--(ws.def("Delete", [](Workspace& w, WorkspaceVariable& x, std::optional<Verbosity> verbosity) {
      Delete(w, Index(1), x.name(), verbosity.has_value() ? *verbosity : WorkspaceVariable{w, )--"
       << verbpos << "});\n}";
    print_method_args(os, *method_ptr, true);
    print_method_desc(os, *method_ptr, arts.varname_group, pass_verbosity);
    os << ");\n\n";

    osptr++;
    if (osptr == oss.end()) osptr = oss.begin();
  }

  for (auto& os : oss) os << "}\n}\n\n";
}

void auto_header(std::ofstream& os, const NameMaps& arts) {
  bool first;

  os << R"--(#ifndef py_auto_interface_h
#define py_auto_interface_h

// Automatic header generated from full Arts setup

#include <agenda_class.h>
#include <pybind11/pybind11.h>
#include <python_interface_value_type.h>
#include <tokval.h>
#include <workspace_ng.h>

namespace Python {
namespace py = pybind11;

using WorkspaceVariablesVariant =
    std::variant<
    )--";

  first = true;
  for (auto& [name, group] : arts.group) {
    if (not first) os << ",\n    ";
    first = false;
    os << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *";
  }

  os << R"--(,
    py::object *>;
  
struct WorkspaceVariablesVariant2String {
)--";

  for (auto& [name, group] : arts.group) {
    os << "constexpr std::string_view operator()(" << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *) const noexcept {return \"" << name << "\";}\n";
  }

  os << R"--(
constexpr std::string_view operator()(Index *) const noexcept {return "Index";}
constexpr std::string_view operator()(Numeric *) const noexcept {return "Numeric";}
constexpr std::string_view operator()(py::object *) const noexcept {return "Python Object";}
};

struct WorkspaceVariable {
  Workspace &ws;
  Index pos;

  WorkspaceVariable(Workspace& ws_, Index pos_) : ws(ws_), pos(pos_) {}

  WorkspaceVariable(Workspace& ws_, Index group_index, const py::object& obj, bool allow_casting);

  [[nodiscard]] const char * name() const;

  [[nodiscard]] Index group() const;

  void pop_workspace_level();

  operator TokVal() const;

  WorkspaceVariable(const WorkspaceVariable&) = default;

  WorkspaceVariable& operator=(const WorkspaceVariable& wv2) {
    pos = wv2.pos;
    return *this;
  }

  [[nodiscard]] Index stack_depth() {return ws.depth(pos);}
  [[nodiscard]] bool is_initialized() const; 
  void initialize_if_not();
  operator WorkspaceVariablesVariant();
  operator WorkspaceVariablesVariant() const;

)--";

  for (auto& [name, group] : arts.group) {
    os << "  operator " << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << "&();\n";
    os << "  operator " << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << "&() const;\n";
  }

  os << R"--(
  // Special cases to allow input to be understood better
  operator Index&();
  operator Index&() const;
  operator Numeric&();
  operator Numeric&() const;
};  // struct WorkspaceVariable
}  // namespace Python

#endif  // py_auto_interface_h
)--";
}

struct TypeVal {
  String type;
  String val;
};

void auto_header_definitions(std::ofstream& os, const NameMaps& arts) {
  os << R"--(#include <python_interface.h>
#include <pybind11/functional.h>

namespace Python {

bool WorkspaceVariable::is_initialized() const { return ws.is_initialized(pos); }

const char * WorkspaceVariable::name() const {
  return ws.wsv_data[pos].Name().c_str();
}

Index WorkspaceVariable::group() const {
  return ws.wsv_data[pos].Group();
}

WorkspaceVariable::WorkspaceVariable(Workspace& ws_, Index group_index, const py::object& obj, bool allow_casting) : ws(ws_), pos(-1) {
  ARTS_USER_ERROR_IF(std::addressof(ws) not_eq std::addressof(ws_), "Cannot use a workspace variable from another workspace")

  if (py::isinstance<WorkspaceVariable>(obj)) {
    pos = obj.cast<WorkspaceVariable *>() -> pos;
  } else {
    ARTS_USER_ERROR_IF(not allow_casting, "Only workspace variable objects allowed")
    ARTS_USER_ERROR_IF(group_index == )--"
     << arts.group.at("Any") << R"--(, "Cannot create type Any")
    std::shared_ptr<void> value_ptr = nullptr;
    switch(group_index) {
)--";

  for (auto& [name, group] : arts.group) {
    char extra = (name == "Numeric" or name == "Index") ? '_' : ' ';
    if (name == "Any") continue;
    os << "      case " << group << ": {auto conv_obj = py::type::of<" << name
       << extra << ">()(";
    if (name == "Agenda") os << "ws, ";
    os << "obj); value_ptr = std::make_shared<" << name << ">(* conv_obj.cast<" << name
       << extra << "*>());} break;\n";
  }

  os << R"--(      default: ARTS_USER_ERROR("Cannot create type")
    }

    // Create the variable and initialize it
    static std::size_t i = 0;
    pos = ws.add_wsv_inplace(WsvRecord(var_string("::anon::", i++).c_str(), "Anonymous agenda variable", group_index));
    ws.push_move(pos, std::move(value_ptr));
  }
}

void WorkspaceVariable::initialize_if_not() {
  if (not ws.is_initialized(pos)) ws.emplace(pos);
  ws[pos];
}

WorkspaceVariable::operator WorkspaceVariablesVariant() {
  initialize_if_not();

  switch (ws.wsv_data[pos].Group()) {
)--";

  for (auto& [name, group] : arts.group) {
    os << "    case " << group << ": return static_cast<" << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>(ws[pos].get());\n";
  }

  os << R"--(  }
  
  ARTS_USER_ERROR("Cannot understand type")
}

WorkspaceVariable::operator WorkspaceVariablesVariant() const {
  ARTS_USER_ERROR_IF(not is_initialized(), "Not initialized: ", ws.wsv_data[pos].Name())

  switch (ws.wsv_data[pos].Group()) {
)--";

  for (auto& [name, group] : arts.group) {
    os << "    case " << group << ": return static_cast<" << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>(ws[pos].get());\n";
  }

  os << R"--(  }
  
  ARTS_USER_ERROR("Cannot understand type")
}

void WorkspaceVariable::pop_workspace_level() { ws.pop(pos); }

  WorkspaceVariable::operator TokVal() const {
  ARTS_USER_ERROR_IF(not is_initialized(), "Must be initialized")
  
  switch (ws.wsv_data[pos].Group()) {
)--";

  for (auto& [name, group] : arts.group) {
    if (name == "Agenda" or name == "ArrayOfAgenda") continue;

    os << "    case " << group << ": return ";
    os << "* static_cast<" << name;
    os << " *>(ws[pos].get())";
    os << ";\n";
  }

  os << R"--(  }
  ARTS_USER_ERROR("Cannot support Agenda and ArrayOfAgenda")
}

Index create_workspace_gin_default_internal(Workspace& ws, const String& key) {
  auto ptr = ws.WsvMap.find(key);
  
  Index pos = ptr == ws.WsvMap.end() ? -1 : ptr -> second;
  )--";

  std::map<String, TypeVal> has;

  for (auto& method : arts.methodname_method) {
    for (std::size_t i = 0; i < method.gin.name.size(); i++) {
      if (method.gin.hasdefs[i]) {
        has[var_string("::", method.name, "::", method.gin.name[i])] =
            TypeVal{method.gin.group[i], method.gin.defs[i]};
      }
    }
  }

  for (auto& [key, items] : has) {
    os << "if (key == \"" << key << "\") {\n";
    os << "    if (pos < 0) pos = ws.add_wsv_inplace(WsvRecord(\""
       << key << R"(", "do not modify", )"
       << arts.group.at(items.type) << "));\n";
    os << "    static_cast<" << items.type
       << " &>(WorkspaceVariable{ws, pos}) = " << items.type << '('
       << items.val << ')';
    os << ";\n";
    os << "  } else ";
  }

  os << "ARTS_USER_ERROR(\"Cannot understand internal key\")\n  return pos;\n}\n\n";

  for (auto& [name, group] : arts.group) {
    os << "WorkspaceVariable::operator " << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << "&() {return *std::get<" << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>(WorkspaceVariablesVariant(*this));}\n";
    os << "WorkspaceVariable::operator " << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << "&() const {return *std::get<" << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>(WorkspaceVariablesVariant(*this));}\n";
  }

  os << R"--(WorkspaceVariable::operator Index&() {return *std::get<Index_ *>(WorkspaceVariablesVariant(*this));}
WorkspaceVariable::operator Index&() const {return *std::get<Index_ *>(WorkspaceVariablesVariant(*this));}
WorkspaceVariable::operator Numeric&() {return *std::get<Numeric_ *>(WorkspaceVariablesVariant(*this));}
WorkspaceVariable::operator Numeric&() const {return *std::get<Numeric_ *>(WorkspaceVariablesVariant(*this));}
}
)--";
}

void workspace_access(std::ofstream& os, const NameMaps& arts) {
  const Index verbpos = arts.varname_group.find("verbosity")->second.artspos;

  for (auto& [name, group] : arts.group) {
    os << R"--(ws.def("_setattr_unchecked_", [](Workspace& w, const char* name, )--"
       << name;
    if (name == "Index" or name == "Numeric") os << '_';
    os << R"--(& val) {
  auto varpos = w.WsvMap.find(name);
  Index i = (varpos == w.WsvMap.end()) ? w.add_wsv_inplace(WsvRecord(name, "User-created variable", )--"
       << group << R"--()) : varpos -> second;
  )--" << name
       << "& x = WorkspaceVariable{w, i};\n";
    if (name == "Agenda")
      os << "  val.set_name(name);\n  val.check(w, *static_cast<Verbosity*>(w["
         << verbpos << "].get()));\n";
    if (name == "ArrayOfAgenda")
      os << "  for (auto& y: val) y.set_name(name);\n  for (auto& y: val) y.check(w, *static_cast<Verbosity*>(w["
         << verbpos << "].get()));\n";
    os << R"--(  x = val;
});

)--";
  }

  os << R"--(wsv.def_property("value",
  py::cpp_function([](const WorkspaceVariable& w) -> WorkspaceVariablesVariant {
    return w; }, py::return_value_policy::reference_internal),
  py::cpp_function([](WorkspaceVariable& w, WorkspaceVariablesVariant v) {
    WorkspaceVariablesVariant wvv{w};
    )--";
  for (auto& [name, group] : arts.group) {
    char extra = (name == "Index" or name == "Numeric") ? '_' : ' ';

    os << "if (std::holds_alternative<" << name << extra << "*>(wvv)) {\n";

    os << "      " << name << "& lh = *std::get<" << name << extra
       << "*>(wvv);\n";

    os << "      if (std::holds_alternative<" << name << extra << "*>(v)) {\n";

    os << "        " << name << "& rh = *std::get<" << name << extra
       << "*>(v);\n";

    if (name == "Agenda")
      os << "        rh.set_name(w.name());\n        rh.check(w.ws, *static_cast<Verbosity*>(w.ws["
         << verbpos << "].get()));\n";
    if (name == "ArrayOfAgenda")
      os << "        for (auto& y: rh) y.set_name(w.name());\n        for (auto& y: rh) y.check(w.ws, *static_cast<Verbosity*>(w.ws["
         << verbpos << "].get()));\n";

    os << "        lh = rh;\n";

    os << "      } else {\n";
    os << "        auto my_obj = py::type::of<" << name << extra << ">()(";
    if (name == "Agenda") os << "w.ws, ";
    os << "* std::get<py::object *>(v));\n";
    os << "        " << name << "& rh = * my_obj.cast<" << name << extra
       << "*>();\n";

    if (name == "Agenda")
      os << "        rh.set_name(w.name());\n        rh.check(w.ws, *static_cast<Verbosity*>(w.ws["
         << verbpos << "].get()));\n";
    if (name == "ArrayOfAgenda")
      os << "        for (auto& y: rh) y.set_name(w.name());\n        for (auto& y: rh) y.check(w.ws, *static_cast<Verbosity*>(w.ws["
         << verbpos << "].get()));\n";

    os << "        lh = rh;\n";

    os << "      }\n";
    os << "    } else ";
  }
  os << R"--({}
  }, py::arg("ws").noconvert(), py::arg("value"))
);

)--";

  for (auto& [name, group] : arts.group) {
    os << "py::implicitly_convertible<WorkspaceVariable, " << name;
    if (name == "Index" or name == "Numeric") os << '_';
    os << ">();\n";
  }
  os << '\n';
}

void make_workspace_pickle(std::ofstream& os, const NameMaps& arts) {
  os << R"--(
ws.def(py::pickle(
  [](Workspace& w) {
    constexpr Index ag_index = )--"
     << arts.group.at("Agenda");
  os << ";\n    constexpr Index aag_index = " << arts.group.at("ArrayOfAgenda");
  os << ";\n    constexpr Index arq_index = "
     << arts.group.at("ArrayOfRetrievalQuantity");
  os << ";\n    constexpr Index cf_index = "
     << arts.group.at("CallbackFunction");
  os << R"--(;
    
    std::vector<py::object> value;
    std::vector<std::string> name;
    std::vector<std::string> group;
    std::vector<std::string> description;

    for (Index i=0; i<w.nelem(); i++) {
      auto& wsv_data = w.wsv_data[i];
      if (w.is_initialized(i) and             // is initialized
          wsv_data.Group() != ag_index and    // is not Agenda (needs workspace variables in fixed positions)
          wsv_data.Group() != aag_index and   // is not ArrayOfAgenda (as above)
          wsv_data.Group() != cf_index and    // is not CallbackFunction (needs code)
          wsv_data.Group() != arq_index and   // is not ArrayOfRetrievalQuantity (always set together with agenda)
          wsv_data.Name().front() != ':') {   // is not generated variable
        name.push_back(wsv_data.Name());
        description.push_back(wsv_data.Description());
        switch (wsv_data.Group()) {)--";
  for (auto& [name, id] : arts.group) {
    const char extra = name == "Index" or name == "Numeric" ? '_' : ' ';
    os << "          case " << id << ": {\n";
    os << "            group.emplace_back(\"" << name << "\");\n";
    os << "            " << name << extra << " val{*static_cast<" << name
       << " *>(w[i].get())};\n";
    os << "            value.emplace_back(py::cast(val));\n";
    os << "          } break;\n";
  }
  os << R"--(          default: { ARTS_USER_ERROR("Bad group") }
        }
      }
    }

    return py::make_tuple(group, name, description, value);
  },
  [](const py::tuple& t) {
    ARTS_USER_ERROR_IF(t.size() != 4, "Invalid state!")

    auto group = t[0].cast<std::vector<std::string>>();
    auto name = t[1].cast<std::vector<std::string>>();
    auto description = t[2].cast<std::vector<std::string>>();
    auto value = t[3].cast<std::vector<py::object>>();
    std::size_t n = group.size();

    ARTS_USER_ERROR_IF(n != name.size() or n != description.size() or n != value.size(), "Invalid state!")

    auto* out = new Workspace{};
    out -> initialize();

    for (std::size_t i=0; i<n; i++) {
      auto wsv_data = out -> WsvMap.find(name[i]);
      const Index wsv_pos =
        (wsv_data == out -> WsvMap.end()) ?
          out -> add_wsv_inplace(WsvRecord(name[i].c_str(), description[i].c_str(), group[i])) :
          wsv_data -> second;

      switch (out -> wsv_data[wsv_pos].Group()) {
)--";
  for (auto& [name, id] : arts.group) {
    const char extra = name == "Index" or name == "Numeric" ? '_' : ' ';
    os << "        case " << id << ": {\n";
    os << "          out -> push_move(wsv_pos, std::make_shared<" << name << ">(* value[i].cast<"
       << name << extra << " *>()));\n";
    os << "          } break;\n";
  }
  os << R"--(        default: { ARTS_USER_ERROR("Bad group") }
      }
    }

    return out;
  }
));
)--";
}

int main(int argc, char** argv) {
  define_wsv_groups();
  Workspace::define_wsv_data();
  Workspace::define_wsv_map();
  define_md_data_raw();
  expand_md_data_raw_to_md_data();
  define_md_map();
  define_agenda_data();
  define_agenda_map();

  ARTS_USER_ERROR_IF(argc not_eq 5,
                     "Usage: EXE NUM_WSM NUM_WSV NUM_CREATE NUM_GENERIC")
  int nwsm = std::stoi(argv[1]);
  int nwsv = std::stoi(argv[2]);
  int nwsc = std::stoi(argv[3]);
  int nwsg = std::stoi(argv[4]);
  ARTS_USER_ERROR_IF(nwsm < 1 or nwsv < 1 or nwsc < 1 or nwsg < 1,
                     "Bad value: ",
                     nwsm,
                     ' ',
                     nwsv,
                     ' ',
                     nwsc,
                     ' ',
                     nwsg)

  const auto artsname = NameMaps();
  {
    std::ofstream py_header("py_auto_interface.h");
    auto_header(py_header, artsname);
  }
  {
    std::ofstream py_header_cc("py_auto_interface.cc");
    auto_header_definitions(py_header_cc, artsname);
  }

  std::ofstream py_workspace("py_auto_workspace.cc");
  includes(py_workspace);
  py_workspace << '\n' << "namespace Python {\n";
  for (int i = 0; i < nwsv; i++) {
    py_workspace << "void py_auto_workspace_wsv_" << i
                 << "(py::class_<Workspace>& ws);\n";
  }
  for (int i = 0; i < nwsg; i++) {
    py_workspace << "void py_auto_workspace_wsg_" << i
                 << "(py::class_<Workspace>& ws);\n";
  }
  for (int i = 0; i < nwsc; i++) {
    py_workspace << "void py_auto_workspace_wsc_" << i
                 << "(py::class_<Workspace>& ws);\n";
  }
  for (int i = 0; i < nwsm; i++) {
    py_workspace << "void py_auto_workspace_wsm_" << i
                 << "(py::class_<Workspace>& ws);\n";
  }
  py_workspace
      << "void py_auto_workspace(py::class_<Workspace>& ws, py::class_<WorkspaceVariable>& wsv) {\n";
  workspace_access(py_workspace, artsname);
  for (int i = 0; i < nwsv; i++) {
    py_workspace << "py_auto_workspace_wsv_" << i << "(ws);\n";
  }
  for (int i = 0; i < nwsg; i++) {
    py_workspace << "py_auto_workspace_wsg_" << i << "(ws);\n";
  }
  for (int i = 0; i < nwsm; i++) {
    py_workspace << "py_auto_workspace_wsm_" << i << "(ws);\n";
  }
  for (int i = 0; i < nwsc; i++) {
    py_workspace << "py_auto_workspace_wsc_" << i << "(ws);\n";
  }
  make_workspace_pickle(py_workspace, artsname);
  py_workspace << "}\n}  // namespace Python\n";

  workspace_variables(nwsv, artsname);
  workspace_method_nongenerics(nwsm, artsname);
  workspace_method_generics(nwsg, artsname);
  workspace_method_create(nwsc, artsname);
}
