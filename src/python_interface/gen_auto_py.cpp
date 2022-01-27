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
  os << "#include <auto_md.h>" << '\n'
     << "#include <m_append.h>" << '\n'
     << "#include <m_basic_types.h>" << '\n'
     << "#include <m_conversion.h>" << '\n'
     << "#include <m_copy.h>" << '\n'
     << "#include <m_delete.h>" << '\n'
     << "#include <m_extract.h>" << '\n'
     << "#include <m_general.h>" << '\n'
     << "#include <m_gridded_fields.h>" << '\n'
     << "#include <m_ignore.h>" << '\n'
     << "#include <m_nc.h>" << '\n'
     << "#include <m_reduce.h>" << '\n'
     << "#include <m_select.h>" << '\n'
     << "#include <m_xml.h>" << '\n'
     << "#include <xml_io.h>" << '\n'
     << '\n'
     << "#include <python_interface.h>" << '\n';
}

void workspace_variables(std::array<std::ofstream, num_split_files>& oss, const NameMaps& arts) {
  auto * osptr = oss.begin();
  for (auto& [name, data] : arts.varname_group) {
    auto& os = *osptr;

    os << "  ws.def_property(\"" << name << "\",\n";
    os << "    [](Workspace& w_) ";
    os << "{ARTS_USER_ERROR_IF (not w_.is_initialized("
       << data.artspos << "), \"Not initialized: " << name
       << "\") return * reinterpret_cast<" << data.varname_group;
    if (data.varname_group == "Index" or data.varname_group == "Numeric") os << '_';
    os << " *>(w_[" << data.artspos << "]);},\n";
    os << "    [](Workspace& w_, " << data.varname_group
       << " val) {if (not w_.is_initialized(" << data.artspos << ")) {w_.push("
       << data.artspos << ", new " << data.varname_group
       << "{});} (* reinterpret_cast<" << data.varname_group << " *>(w_["
       << data.artspos << "])) = std::move(val);},\n";
    os << "    R\"-VARNAME_DESC-(" << data.varname_desc << ")-VARNAME_DESC-\""
       << ");";
    os << '\n' << '\n';

    osptr++;
    if (osptr == oss.end()) osptr = oss.begin();
  }
}

void workspace_methods(std::array<std::ofstream, num_split_files>& oss, const NameMaps& arts) {
  auto * osptr = oss.begin();

  for (auto& method : arts.methodname_method) {
    auto& os = *osptr;

    if (method.agenda_method) continue;  // FIXME: Should be fixed???

    // Skip create methods since these must be called via Var
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

    const bool has_many = [&](){
      long i=0;
      for (auto& a: arts.methodname_method)
        i += (a.name == method.name);
      return i > 1;
    }();

    const bool is_last = [&](){
      auto * last = &method;
      for (auto& a : arts.methodname_method) {
        if (a.name == method.name) last = &a;
      }
      return last == &method;
    }();

    // Arguments from python side
    os << "ws.def(\"" << method.name << "\",\n[](Workspace& w_ [[maybe_unused]], ";
    for (const auto & i : method.out.varname) {
      auto& group = arts.varname_group.at(i).varname_group;
      os << group;
      if (group == "Index" or group == "Numeric") os << "_";
      os << "* " << i << ", ";
    }
    for (std::size_t i = 0; i < method.gout.name.size(); i++) {
      auto& group = method.gout.group[i];
      os << group;
      if (group == "Index" or group == "Numeric") os << "_";
      os << "* " << method.gout.name[i] << ", ";
    }
    for (const auto & i: method.in.varname) {
      if (std::none_of(method.out.varname.cbegin(),
                       method.out.varname.cend(),
                       [in = i](const auto& out) { return in == out; })) {
        auto& group = arts.varname_group.at(i).varname_group;
        os << group;
        if (group == "Index" or group == "Numeric") os << "_";
        os << "* " << i << ", ";
      }
    }
    for (std::size_t i = 0; i < method.gin.name.size(); i++) {
      auto& group = method.gin.group[i];
      os << group;
      if (group == "Index" or group == "Numeric") os << "_";
      os << "* " << method.gin.name[i] << ", ";
    }
    if (pass_verbosity) os << "Verbosity* verbosity, ";
    os << "py::kwargs& kwargs) {\n";
    bool has_any = false;

    // Create static defaults
    for (std::size_t i = 0; i < method.gin.name.size(); i++) {
      if (method.gin.hasdefs[i]) {
        os << "static " << method.gin.group[i] << " " << method.gin.name[i] << "_{";
        if (method.gin.defs[i] not_eq "{}") os << method.gin.defs[i];
        os << "};\n";
      }
    }

    // Pure outputs need a stack
    for (const auto & i : method.out.varname) {
      if (std::none_of(method.in.varname.cbegin(),
                       method.in.varname.cend(),
                       [in = i](const auto& out) { return in == out; })) {
        auto& group = arts.varname_group.at(i).varname_group;
        auto& pos = arts.varname_group.at(i).artspos;
        os << "if (not w_.is_initialized(" << pos << ")) w_.push("<<pos<<", new "<<group<<"{});\n";
      }
    }

    // Arguments from Arts side
    os << method.name << '(';
    if (pass_workspace) {
      os << "w_";
      has_any = true;
    } else os <<'\n';

    // Outputs
    for (const auto & i : method.out.varname) {
      if (has_any) os << ",\n";
      has_any = true;

      os << "select<";
      auto& group = arts.varname_group.at(i).varname_group;
      os << group;
      os << ">(w_, " << arts.varname_group.at(i).artspos << ", " << i
         << ", kwargs, \"" << i << "\")";
    }

    // Generic Outputs
    for (std::size_t i = 0; i < method.gout.name.size(); i++) {
      if (has_any) os << ",\n";
      has_any = true;

      os << "select<";
      os << method.gout.group[i];
      os << ">(" << method.gout.name[i] << ", kwargs, \"" << method.gout.name[i]
         << "\")";
    }

    // False Generic Output Names
    if (method.pass_wsv_names) {
      for (auto& name : method.gout.name) {
        if (has_any) os << ",\n";
        has_any = true;
        os << "\"" << name << "\"";
      }
    }

    // Inputs (that are not outputs)
    for (auto& i: method.in.varname) {
      if (std::none_of(method.out.varname.cbegin(),
                       method.out.varname.cend(),
                       [in = i](const auto& out) { return in == out; })) {
        if (has_any) os << ",\n";
        has_any = true;

        os << "select<";
        auto& group = arts.varname_group.at(i).varname_group;
        os << group;
        os << ">(w_, " << arts.varname_group.at(i).artspos << ", " << i
          << ", kwargs, \"" << i << "\")";
      }
    }

    // Generic Inputs
    for (std::size_t i = 0; i < method.gin.name.size(); i++) {
      if (has_any) os << ",\n";
      has_any = true;
      os << "select<";
      os << method.gin.group[i];
      os << ">(";
      if (method.gin.hasdefs[i]) os << method.gin.name[i] << "_, ";
      os << method.gin.name[i] << ", kwargs, \"" << method.gin.name[i]
         << "\")";
    }

    // False Generic Input Names
    if (method.pass_wsv_names) {
      for (auto& name : method.gin.name) {
        if (has_any) os << ",\n";
        has_any = true;
        os << "\"" << name << "\"";
      }
    }

    // Verbosity?
    if (pass_verbosity) {
      if (has_any) os << ",\n";
      os <<"select<Verbosity>(w_, " << arts.varname_group.at("verbosity").artspos << ", verbosity, kwargs, \"verbosity\")";
    }
    os << ");\n}";

    for (const auto& i : method.out.varname) {
      os << ','  << '\n' << "py::arg_v(\"" << i << "\", nullptr, \"Workspace::" << i << "\")";
    }
    for (const auto& i : method.gout.name) {
      os << ','  << '\n' << "py::arg(\"" << i << "\") = nullptr";
    }
    for (const auto& i : method.in.varname) {
      if (std::none_of(method.out.varname.cbegin(),
                       method.out.varname.cend(),
                       [in = i](const auto& out) { return in == out; })) {
        os << ','  << '\n' << "py::arg_v(\"" << i << "\", nullptr, \"Workspace::" << i << "\")";
      }
    }
    for (size_t i = 0; i<method.gin.name.size(); i++) {
      if (method.gin.hasdefs[i]) {
        if (method.gin.group[i] == "String")
          os << ','  << '\n' << "py::arg_v(\"" << method.gin.hasdefs[i] << R"(", nullptr, "\)" << method.gin.defs[i] << R"("\"")" << ')';
        else 
          os << ','  << '\n' << "py::arg_v(\"" << method.gin.hasdefs[i] << "\", nullptr, \"" << method.gin.defs[i] << "\")";
      } else {
        os << ','  << '\n' << "py::arg(\"" << method.gin.name[i] << "\") = nullptr";
      }
    }
    if (pass_verbosity) os << ','  << '\n' << R"--(py::arg_v("verbosity", nullptr, "Workspace::verbosity"))--";

    // Put comment at the end
    if (not has_many or is_last) {
      os << ",\npy::doc(\nR\"-METHODS_DESC-(\n" << method.desc;
      os << "\n\nAuthor";
      if (method.authors.size() > 1) os << "s";
      os << ": ";
      for (auto& author: method.authors) os << author << ", ";
      os << "\n\")-METHODS_DESC-\")";
    }

    os << ')' << ';' << '\n' << '\n';

    if (not has_many or is_last) {
      osptr++;
      if (osptr == oss.end()) osptr = oss.begin();
    }
  }
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

  std::ofstream py_workspace("py_auto_workspace.cc");
  std::array<std::ofstream, num_split_files> py_workspaces;
  for (Index i=0; i<num_split_files; i++) py_workspaces[i] = std::ofstream(var_string("py_auto_workspace_split_", i, ".cc"));

  includes(py_workspace);
  for (Index i=0; i<num_split_files; i++) includes(py_workspaces[i]);

  const auto artsname = NameMaps();

  py_workspace << '\n' << "namespace Python {\n";
  for (Index i=0; i<num_split_files; i++) py_workspace << "void py_auto_workspace_" << i << "(py::class_<Workspace>& ws);\n";
  py_workspace << "void py_auto_workspace(py::class_<Workspace>& ws) {\n";
  for (Index i=0; i<num_split_files; i++) py_workspaces[i] << '\n' << "namespace Python {\nvoid py_auto_workspace_" << i << "(py::class_<Workspace>& ws) {\n";
  workspace_variables(py_workspaces, artsname);
  workspace_methods(py_workspaces, artsname);
  for (Index i=0; i<num_split_files; i++) py_workspace << "py_auto_workspace_" << i << "(ws);\n";
  py_workspace << "}\n}  // namespace Python\n";
  for (Index i=0; i<num_split_files; i++) py_workspaces[i] << "}\n}  // namespace Python\n";
}
