#include <auto_md.h>
#include <global_data.h>

#include <algorithm>
#include <fstream>
#include <ostream>
#include <utility>

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
  for (auto& x : Workspace::wsv_data) desc[x.Name()] = x.Description();
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
      defaults[i] = "{}";
    } else {
      defaults[i] = vardefaults[i];
    }

    if (defaults[i] == "") defaults[i] = "{}";

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

// To speed up compilation, we split the workspace into N units
constexpr Index num_split_files = 32;

void includes(std::ofstream& os) {
  os << "#include <py_auto_interface.h>" << '\n';
}

void workspace_variables(std::array<std::ofstream, num_split_files>& oss,
                         const NameMaps& arts) {
  auto* osptr = oss.begin();
  for (auto& [name, data] : arts.varname_group) {
    auto& os = *osptr;

    os << "  ws.def_property(\"" << name << "\",\n";
    os << "    py::cpp_function([](Workspace& w) -> WorkspaceVariable ";
    os << "{return WorkspaceVariable{w, " << data.artspos
       << "};}, py::return_value_policy::reference_internal),\n";
    os << "    [](Workspace& w, " << data.varname_group
       << " val) {"<<data.varname_group<<"& v_ = WorkspaceVariable{w, " << data.artspos << "}; v_ = std::move(val);},\n";
    os << "    R\"-VARNAME_DESC-(\n"
       << data.varname_desc << ")-VARNAME_DESC-\""
       << ");";
    os << '\n' << '\n';

    osptr++;
    if (osptr == oss.end()) osptr = oss.begin();
  }
}

void workspace_method_nongenerics(std::array<std::ofstream, num_split_files>& oss,
                       const NameMaps& arts) {
  auto* osptr = oss.begin();

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
    }()) continue;  // Skip 'generic' functions

    // Arguments from python side
    os << "ws.def(\"" << method.name
       << "\",[](Workspace& w_ [[maybe_unused]]";
    for (const auto& i : method.out.varname) {
      auto& group = arts.varname_group.at(i).varname_group;
      os << ",\nstd::optional<std::variant<";
      if (std::any_of(method.in.varname.cbegin(),
                      method.in.varname.cend(),
                      [out = i](const auto& in) { return in == out; })) os << "const ";
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
        os << ",\nstd::optional<std::variant<const WorkspaceVariable *, " << group;
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
    if (pass_verbosity) os << ",\nstd::optional<std::variant<const WorkspaceVariable *, Verbosity *>> verbosity";
    os << ") {\n";
    bool has_any = false;

    // Create static defaults
    for (std::size_t i = 0; i < method.gin.name.size(); i++) {
      if (method.gin.hasdefs[i]) {
        os << "static const " << method.gin.group[i] << " " << method.gin.name[i]
           << "_{";
        if (method.gin.defs[i] not_eq "{}") os << method.gin.defs[i];
        os << "};\n";
      }
    }

    // Outputs
    Index counter = 0;
    for (const auto& i : method.out.varname) {
      auto& group = arts.varname_group.at(i).varname_group;
      os << group << "& arg" << counter++ << "_ = select_";
      if (std::any_of(method.in.varname.cbegin(),
                      method.in.varname.cend(),
                      [out = i](const auto& in) { return in == out; })) os << "in";
      os << "out<";
      os << group;
      os << ">(WorkspaceVariable{w_, " << arts.varname_group.at(i).artspos << "}, " << i << ");\n";
    }

    // Generic Outputs
    for (std::size_t i = 0; i < method.gout.name.size(); i++) {
      os << method.gout.group[i] << "& arg" << counter++ << "_ = select_gout<";
      os << method.gout.group[i];
      os << ">(" << method.gout.name[i] << ");\n";
    }

    // False Generic Output Names
    if (method.pass_wsv_names) {
      for (auto& name : method.gout.name) {
        os << "const String arg" << counter++ << "_{\"" << name << "\"};\n";
      }
    }

    // Inputs (that are not outputs)
    for (auto& i : method.in.varname) {
      if (std::none_of(method.out.varname.cbegin(),
                       method.out.varname.cend(),
                       [in = i](const auto& out) { return in == out; })) {
        auto& group = arts.varname_group.at(i).varname_group;
        os << "const "<<group<<"& arg" << counter++ << "_ = select_in<";
        os << group;
        os << ">(WorkspaceVariable{w_, " << arts.varname_group.at(i).artspos << "}, " << i
           << ");\n";
      }
    }

    // Generic Inputs
    for (std::size_t i = 0; i < method.gin.name.size(); i++) {
      os << "const " << method.gin.group[i] << "& arg" << counter++
         << "_ = select_gin<" << method.gin.group[i] << ">(";
      if (method.gin.hasdefs[i]) os << method.gin.name[i] << "_, ";
      os << method.gin.name[i] << ");\n";
    }

    // False Generic Input Names
    if (method.pass_wsv_names) {
      for (auto& name : method.gin.name) {
        os << "const String arg" << counter++ << "_{\"" << name << "\"};\n";
      }
    }

    // Verbosity?
    if (pass_verbosity) {
      os << "const Verbosity& arg" << counter++ << "_ = select_in<Verbosity>(WorkspaceVariable{w_, "
         << arts.varname_group.at("verbosity").artspos << "}, verbosity);\n";
    }

    // Arguments from Arts side
    has_any = false;
    os << method.name << '(';
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
    for (const auto& i : method.out.varname) {
      os << ',' << '\n'
         << "py::arg_v(\"" << i << "\", std::nullopt, \"Workspace::" << i << "\")";
    }
    for (const auto& i : method.gout.name) {
      os << ',' << '\n' << "py::arg(\"" << i << "\")";
    }
    for (const auto& i : method.in.varname) {
      if (std::none_of(method.out.varname.cbegin(),
                       method.out.varname.cend(),
                       [in = i](const auto& out) { return in == out; })) {
        os << ',' << '\n'
           << "py::arg_v(\"" << i << "\", std::nullopt, \"Workspace::" << i << "\")";
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
        os << ',' << '\n' << "py::arg(\"" << method.gin.name[i] << "\")";
      }
    }
    if (pass_verbosity)
      os << ',' << '\n'
         << R"--(py::arg_v("verbosity", std::nullopt, "Workspace::verbosity"))--";

    // Put description at the end
    os << ",\npy::doc(\nR\"-METHODS_DESC-(\n" << method.desc;
    os << "\n\nAuthor";
    if (method.authors.size() > 1) os << "s";
    os << ": ";
    for (auto& author : method.authors) os << author << ", ";
    os << "\n)-METHODS_DESC-\")";

    os << ')' << ';' << '\n' << '\n';

    osptr++;
    if (osptr == oss.end()) osptr = oss.begin();
  }
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

void workspace_method_generics(std::array<std::ofstream, num_split_files>& oss,
                       const NameMaps& arts) {
  auto* osptr = oss.begin();

  for (std::size_t imethod=0; imethod<arts.methodname_method.size(); imethod++) {
    auto& method = arts.methodname_method[imethod];
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
                               [out = meth.out.varname[i]](const auto& in) { return in == out; });
          arg.out = true;
          arg.gen = false;
          arg.def = false;
          arg.wsv = true;
          arg.wsvpos = meth.out.varpos[i];
          arg.name = meth.out.varname[i];
          arg.def_str = "";
          arg.types = ArrayOfString{arts.varname_group.at(arg.name).varname_group};
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
            arg.wsvpos = meth.in.varpos[i];
            arg.name = meth.in.varname[i];
            arg.def_str = "";
            arg.types = ArrayOfString{arts.varname_group.at(arg.name).varname_group};
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
      first=false;
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

    bool has_any = false;

    // Arguments from python side
    os << "ws.def(\"" << method.name
       << "\",[](Workspace& w_ [[maybe_unused]]";
    for (auto& arg: arg_help) {
      os << ",\n";
      if (arg.wsv or arg.def) os << "std::optional<";
      os << "std::variant<";
      if (arg.in) os << "const ";
      os << "WorkspaceVariable *, ";
      if (arg.types.size() == 1) {
        os << arg.types.front();
        if (arg.types.front() == "Index" or arg.types.front() == "Numeric") os << "_";
        os << " *";
      } else {
        ArrayOfString copy_types = arg.types;
        std::sort(copy_types.begin(), copy_types.end());
        copy_types.erase(std::unique(copy_types.begin(), copy_types.end()), copy_types.end());

        has_any = false;
        for (auto& t: copy_types) {
          if (has_any) os << ", ";
          has_any = true;
          os << t;
        if (t == "Index" or t == "Numeric") os << "_";
          os << " *";
        }
      }
      os << '>';

      if (not arg.gen or arg.def) os << '>';
      os << ' ' << arg.name;
    }
    if (pass_verbosity) os << ",\nstd::optional<std::variant<const WorkspaceVariable *, Verbosity *>> verbosity";
    os << ") {\n";

    // Create static defaults
    for (std::size_t i = 0; i < method.gin.name.size(); i++) {
      if (method.gin.hasdefs[i]) {
        os << "static const " << method.gin.group[i] << " " << method.gin.name[i]
           << "_{";
        if (method.gin.defs[i] not_eq "{}") os << method.gin.defs[i];
        os << "};\n";
      }
    }

    // Variables
    bool start_gout=false, done_gout_names=false;
    Index counter = 0;
    for (auto& arg: arg_help) {

      // Pass all the names in order
      if (method.pass_wsv_names) {
        if (arg.gen and arg.out) {
          start_gout = true;
        } else if (start_gout and not done_gout_names) {
          for (auto& arg2: arg_help) {
            if (arg2.gen and arg2.out)
              os << "static const String arg" << counter++ << "_{\"" << arg2.name << "\"};\n";
          }
          done_gout_names = true;
        }
      }

      if (arg.types.size() == 1) {
        if (arg.in and not arg.out) os << "const ";
        os << arg.types.front() << "& ";
      } else {
        os << "WorkspaceVariablesVariant wvv_";
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
        os << "wvv<WorkspaceVariablesVariant>(" << arg.name << ");";
      }
      os << '\n';
    }
    if (method.pass_wsv_names) {
      for (auto& arg : arg_help) {
        if (arg.gen and arg.in)
          os << "static const String arg" << counter++ << "_{\"" << arg.name
             << "\"};\n";
      }
    }
    if (pass_verbosity) {
      os << "const Verbosity& arg" << counter++ << "_ = select_in<Verbosity>(WorkspaceVariable{w_, "
         << arts.varname_group.at("verbosity").artspos << "}, verbosity);\n";
    }
    const Index total_count = counter;

    // Find the right function
    for (std::size_t i=0; i< n_func; i++) {
      os << "if (";
      counter = 0;
      has_any = false;

      start_gout = false, done_gout_names=false;
      for (auto& arg: arg_help) {
        if(method.pass_wsv_names) {
          if (arg.gen and arg.out) {
            start_gout = true;
          } else if (start_gout and not done_gout_names) {
            for (auto& arg2: arg_help) {
              if (arg2.gen and arg2.out)
                counter++;
            }
            done_gout_names = true;
          }
        }

        if (arg.types.size() > 1) {
          if (has_any) os << " and ";
          has_any = true;
          os << "std::holds_alternative<" << arg.types[i];
          if(arg.types[i] == "Numeric" or arg.types[i] == "Index") os << '_';
          os << " *>(wvv_arg" << counter << "_)";
        }
        counter++;
      }
      os << ") {\n";

      // Interpret the right variables
      counter = 0;
      start_gout = false, done_gout_names=false;
      for (auto& arg: arg_help) {
        if(method.pass_wsv_names) {
          if (arg.gen and arg.out) {
            start_gout = true;
          } else if (start_gout and not done_gout_names) {
            for (auto& arg2: arg_help) {
              if (arg2.gen and arg2.out)
                counter++;
            }
            done_gout_names = true;
          }
        }

        if (arg.types.size() > 1) {
          if (not arg.out) os << "const ";
          os << arg.types[i];
          os << "& arg" << counter << "_ = *std::get<" << arg.types[i];
          if(arg.types[i] == "Numeric" or arg.types[i] == "Index") os << '_';
          os << " *>(wvv_arg" << counter << "_);\n";
        }
        counter++;
      }
      
        // Arguments from Arts side
    has_any = false;
    os << method.name << '(';
    if (pass_workspace) {
      os << "w_";
      has_any = true;
    }
    for (Index current_count = 0; current_count < total_count; current_count++) {
      if (has_any) os << ", ";
      has_any = true;
      os << "arg" << current_count << "_";
    }
    os << ");\n";
      
      os << "} else ";
    }

    os << "ARTS_USER_ERROR(\"Generic Arts function called but without exact match of input parameter types\")}";

    // Name the paramters and show their defaults
    for (const auto& i : method.out.varname) {
      os << ',' << '\n'
         << "py::arg_v(\"" << i << "\", std::nullopt, \"Workspace::" << i << "\")";
    }
    for (const auto& i : method.gout.name) {
      os << ',' << '\n' << "py::arg(\"" << i << "\")";
    }
    for (const auto& i : method.in.varname) {
      if (std::none_of(method.out.varname.cbegin(),
                       method.out.varname.cend(),
                       [in = i](const auto& out) { return in == out; })) {
        os << ',' << '\n'
           << "py::arg_v(\"" << i << "\", std::nullopt, \"Workspace::" << i << "\")";
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
        os << ',' << '\n' << "py::arg(\"" << method.gin.name[i] << "\")";
      }
    }
    if (pass_verbosity)
      os << ',' << '\n'
         << R"--(py::arg_v("verbosity", std::nullopt, "Workspace::verbosity"))--";

    // Put description at the end
    os << ",\npy::doc(\nR\"-METHODS_DESC-(\n" << method.desc;
    os << "\n\nAuthor";
    if (method.authors.size() > 1) os << "s";
    os << ": ";
    for (auto& author : method.authors) os << author << ", ";
    os << "\n)-METHODS_DESC-\")";

    os << ')' << ';' << '\n' << '\n';

    osptr++;
    if (osptr == oss.end()) osptr = oss.begin();
  }
}

void auto_header(std::ofstream& os, const NameMaps& arts) {
  bool first;

  os << R"--(// Automatic header generated from full Arts setup

#include <auto_md.h>
#include <m_append.h>
#include <m_basic_types.h>
#include <m_conversion.h>
#include <m_copy.h>
#include <m_delete.h>
#include <m_extract.h>
#include <m_general.h>
#include <m_gridded_fields.h>
#include <m_ignore.h>
#include <m_nc.h>
#include <m_reduce.h>
#include <m_select.h>
#include <m_xml.h>
#include <python_interface.h>
#include <supergeneric.h>
#include <xml_io.h>
namespace Python {
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

struct WorkspaceVariable {
  Workspace &ws;
  Index pos;

  WorkspaceVariable(Workspace& ws_, Index pos_) : ws(ws_), pos(pos_) {}

  WorkspaceVariable(Workspace& ws_, Index group_index, py::object obj);

  String name() const;

  WorkspaceVariable& operator=(WorkspaceVariable& wv2) {
    pos = wv2.pos;
    return *this;
  }

  bool is_initialized() const; 
  void initialize_if_not();
  operator WorkspaceVariablesVariant();
  operator WorkspaceVariablesVariant() const;
  WorkspaceVariable &operator=(WorkspaceVariablesVariant x);
)--";

  for (auto& [name, group] : arts.group) {
    os << "  operator " << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os<< "&() {return *std::get<" << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>(WorkspaceVariablesVariant(*this));}\n";
    os << "  operator " << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << "&() const {return *std::get<" << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>(WorkspaceVariablesVariant(*this));}\n";
  }

os << R"--(
  // Special cases to allow input to be understood better
  operator Index&() {return *std::get<Index_ *>(WorkspaceVariablesVariant(*this));}
  operator Index&() const {return *std::get<Index_ *>(WorkspaceVariablesVariant(*this));}
  operator Numeric&() {return *std::get<Numeric_ *>(WorkspaceVariablesVariant(*this));}
  operator Numeric&() const {return *std::get<Numeric_ *>(WorkspaceVariablesVariant(*this));}
)--";

  os << R"--(};  // struct WorkspaceVariable
}  // namespace Python
)--";
}

void auto_header_definitions(std::ofstream& os, const NameMaps& arts) {
  os << R"--(
bool WorkspaceVariable::is_initialized() const { return ws.is_initialized(pos); }

String WorkspaceVariable::name() const {
  return ws.wsv_data[pos].Name();
}

WorkspaceVariable::WorkspaceVariable(Workspace& ws_, Index group_index, py::object obj) : ws(ws_), pos(-1) {
  static std::size_t i = 0;

  try {
    *this = *obj.cast<WorkspaceVariable *>();
    return;
  } catch (...) {
    // Do nothing
  }

  ARTS_USER_ERROR_IF(group_index == )--" << arts.group.at("Any") <<  R"--(, "Cannot create type Any")
  void * value_ptr = nullptr;
  )--";

  for (auto& [name, group] : arts.group) {
    os << "if (group_index == " << group << ") {\n";
    os << "   value_ptr = " << "obj.cast<" << name;
    if (name == "Index" or name == "Numeric") os << '_';
    os << " *>();\n" << "  } else ";
  }

  os << R"--(ARTS_USER_ERROR("Cannot create type")

  // Create the variable and initialize it
  pos = ws.add_wsv_inplace(WsvRecord(var_string("_anon", i++).c_str(), "Created by pybind11 API", group_index));
  ws.push(pos, value_ptr);
}

void WorkspaceVariable::initialize_if_not() {
  if (not ws.is_initialized(pos)) {
    switch (Workspace::wsv_data[pos].Group()) {
)--";

  for (auto& [name, group] : arts.group)
    os << "      case " << group << ":  ws.push(pos, new " << name
       << "); break;\n";

  os << R"--(
      default: ARTS_USER_ERROR("Unknown variable group")
    }
  }
}

WorkspaceVariable::operator WorkspaceVariablesVariant() {
  initialize_if_not();

  switch (Workspace::wsv_data[pos].Group()) {
)--";

  for (auto& [name, group] : arts.group) {
    os << "    case " << group << ": return reinterpret_cast<" << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>(ws[pos]);\n";
  }

  os << R"--(  }
  
  ARTS_USER_ERROR("Cannot understand type")
}

WorkspaceVariable::operator WorkspaceVariablesVariant() const {
  ARTS_USER_ERROR_IF(not is_initialized(), "Not initialized: ", Workspace::wsv_data[pos].Name())

  switch (Workspace::wsv_data[pos].Group()) {
)--";

  for (auto& [name, group] : arts.group) {
    os << "    case " << group << ": return reinterpret_cast<" << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>(ws[pos]);\n";
  }

  os << R"--(  }
  
  ARTS_USER_ERROR("Cannot understand type")
}

WorkspaceVariable &WorkspaceVariable::operator=(WorkspaceVariablesVariant x) {
  initialize_if_not();
  switch (Workspace::wsv_data[pos].Group()) {
)--";
  for (auto& [name, group] : arts.group) {
    os << "    case " << group
       << ":\n"
          "      if (std::holds_alternative<py::object *>(x)) {\n"
          "        * reinterpret_cast<"
       << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>(ws[pos]) = * std::get<py::object *>(x) -> cast<" << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>();\n";

    os << "      } else {\n"
          "        ARTS_USER_ERROR_IF(not std::holds_alternative<"
       << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>(x), \"Cannot cast between internal classes\")\n"
          "        * reinterpret_cast<"
       << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>(ws[pos]) = * std::get<" << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << " *>(x);\n    } break;\n";
  }
  os << R"--(  default: ARTS_USER_ERROR("Unknown variable group: ", Workspace::wsv_data[pos].Name())
  }
  
  return *this;
}
)--";
}

void workspace_access(std::ofstream& os, const NameMaps& arts) {
  os << R"--(
ws.def("__setattr__", [](Workspace& w, const char * name, WorkspaceVariablesVariant x) {
  auto varpos = w.WsvMap.find(name);

  bool newname = varpos == w.WsvMap.end();
  Index i=-1;
  if (newname) {
    ARTS_USER_ERROR_IF(auto _sv=std::string_view(name); _sv.size() == 0 or _sv.front() == '_' or _sv.back() == '_',
    "Cannot automatically generate names starting or ending with \'_\', it is reserved for internal use.\n"
    "If you really need it, use create_variable at your own risk")

    )--";
  bool first = true;
  for (auto& [name, group] : arts.group) {
    if (not first) os << "    else ";
    first = false;
    os << "if (std::holds_alternative<" << name;
    if (name == "Index" or name == "Numeric") os << "_";
    os << R"--( *>(x)) i = w.add_wsv_inplace(WsvRecord(name, "User-generated value", )--"
       << group << "));\n";
  }

  os << R"--(    else
      ARTS_USER_ERROR("Cannot recognize workspace variable type\n\n"
      "You cannot initiate a workspace variable from a pure python object (use create_variable if unsure how to proceed)")
  } else i = varpos->second;

  WorkspaceVariable var{w, i};
  var.initialize_if_not();
  var = x;
});

)--";


  for (auto& [name, group] : arts.group) {
    os << "py::implicitly_convertible<const WorkspaceVariable, "<< name;
    if (name == "Index" or name == "Numeric") os << "_";
    os <<">();\n";
  }
  os << '\n';
}

struct TypeVal {
  String type;
  String val;
};

void internal_defaults(std::ofstream& os, const NameMaps& arts) {
  std::map<String, TypeVal> has;

  for (auto& method: arts.methodname_method) {
    for (std::size_t i=0; i<method.gin.name.size(); i++) {
      if (method.gin.hasdefs[i]) {
        String gin_key = "GeneratedInputDefault_" + method.name + "_" + method.gin.name[i] + "_";
        has[gin_key] = TypeVal{method.gin.group[i], method.gin.defs[i]};
      }
    }
  }

  os << R"--(
Index create_workspace_gin_default_internal(Workspace& ws, String& key) {
  )--";
  
  for (auto& [key, items]: has) {
    os << "if (key == \"" << key << "\") {\n";
    os << "    WorkspaceVariable wv{ws, ws.add_wsv_inplace(WsvRecord(\""<<key<<"\", \"Created by pybind11 API\", " << arts.group.at(items.type) << "))};\n";
    os << "    wv.initialize_if_not();\n";
    os << "    " << items.type << "& val = wv;\n";
    os << "    val = "<<items.type<<"{";
    if ("{}" not_eq items.val) os << items.val;
    os <<"};\n";
    os << "    return wv.pos;\n";
    os << "  } else ";
  }

  os << "ARTS_USER_ERROR(\"Cannot understand internal key\")\n  return -1;\n}\n\n";
}

int main() {
  define_wsv_group_names();
  Workspace::define_wsv_data();
  Workspace::define_wsv_map();
  define_md_data_raw();
  expand_md_data_raw_to_md_data();
  define_md_map();
  define_agenda_data();
  define_agenda_map();

  const auto artsname = NameMaps();

  std::ofstream py_header("py_auto_interface.h");
  auto_header(py_header, artsname);

  std::ofstream py_workspace("py_auto_workspace.cc");
  std::array<std::ofstream, num_split_files> py_workspaces;
  for (Index i = 0; i < num_split_files; i++)
    py_workspaces[i] =
        std::ofstream(var_string("py_auto_workspace_split_", i, ".cc"));

  includes(py_workspace);
  for (Index i = 0; i < num_split_files; i++) includes(py_workspaces[i]);

  py_workspace << '\n' << "namespace Python {\n";
  auto_header_definitions(py_workspace, artsname);
  internal_defaults(py_workspace, artsname);
  for (Index i = 0; i < num_split_files; i++) {
    py_workspace << "void py_auto_workspace_" << i
                 << "(py::class_<Workspace>& ws);\n";
  }
  py_workspace << "void py_auto_workspace(py::class_<Workspace>& ws) {\n";
  for (Index i = 0; i < num_split_files; i++) {
    py_workspaces[i] << '\n'
                     << "namespace Python {\nvoid py_auto_workspace_" << i
                     << "(py::class_<Workspace>& ws) {\n";
  }
  workspace_variables(py_workspaces, artsname);
  workspace_method_nongenerics(py_workspaces, artsname);
  workspace_method_generics(py_workspaces, artsname);
  for (Index i = 0; i < num_split_files; i++) {
    py_workspace << "py_auto_workspace_" << i << "(ws);\n";
  }
  workspace_access(py_workspace, artsname);
  py_workspace << "}\n}  // namespace Python\n";
  for (Index i = 0; i < num_split_files; i++) {
    py_workspaces[i] << "}\n}  // namespace Python\n";
  }
}
