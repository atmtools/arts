#include <auto_md.h>
#include <global_data.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

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
    gin_group.push_back({x.GInType().cbegin(), x.GInType().cend()});
  std::vector<std::vector<std::string>> gin_names;
  for (auto& x : global_data::md_data)
    gin_names.push_back({x.GIn().cbegin(), x.GIn().cend()});
  std::vector<std::vector<std::string>> gin_defaults;
  for (auto& x : global_data::md_data)
    gin_defaults.push_back({x.GInDefault().cbegin(), x.GInDefault().cend()});
  std::vector<std::vector<std::string>> gin_desc;
  for (auto& x : global_data::md_data)
    gin_desc.push_back(
        {x.GInDescription().cbegin(), x.GInDescription().cend()});
  std::vector<std::vector<std::size_t>> gout_group;
  for (auto& x : global_data::md_data)
    gout_group.push_back({x.GOutType().cbegin(), x.GOutType().cend()});
  std::vector<std::vector<std::string>> gout_names;
  for (auto& x : global_data::md_data)
    gout_names.push_back({x.GOut().cbegin(), x.GOut().cend()});
  std::vector<std::vector<std::string>> gout_desc;
  for (auto& x : global_data::md_data)
    gout_desc.push_back(
        {x.GOutDescription().cbegin(), x.GOutDescription().cend()});
  std::vector<std::vector<std::size_t>> in_wspace;
  for (auto& x : global_data::md_data)
    in_wspace.push_back({x.In().cbegin(), x.In().cend()});
  std::vector<std::vector<std::size_t>> out_wspace;
  for (auto& x : global_data::md_data)
    out_wspace.push_back({x.Out().cbegin(), x.Out().cend()});
  std::vector<std::string> desc;
  for (auto& x : global_data::md_data) desc.push_back(x.Description());
  std::vector<std::vector<std::string>> authors;
  for (auto& x : global_data::md_data)
    authors.push_back({x.Authors().cbegin(), x.Authors().cend()});

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
    inoutvarpos.push_back({x.InOut().cbegin(), x.InOut().cend()});
  std::vector<std::vector<std::size_t>> invarpos;
  for (auto& x : global_data::md_data)
    invarpos.push_back({x.InOnly().cbegin(), x.InOnly().cend()});
  std::vector<std::vector<std::size_t>> outvarpos;
  for (auto& x : global_data::md_data)
    outvarpos.push_back({x.Out().cbegin(), x.Out().cend()});

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

int main() {
  define_wsv_group_names();
  Workspace::define_wsv_data();
  Workspace::define_wsv_map();
  define_md_data_raw();
  expand_md_data_raw_to_md_data();
  define_md_map();
  define_agenda_data();
  define_agenda_map();
  define_species_data();
  define_species_map();

  const auto artsname = NameMaps();

  std::cout << "#ifndef autoarts_h\n"
            << "#define autoarts_h\n"
            << '\n'
            << "#include <auto_md.h>" << '\n'
            << "#include <arts.h>" << '\n'
            << "#include <global_data.h>" << '\n'
            << "#include <m_basic_types.h>" << '\n'
            << "#include <m_general.h>" << '\n'
            << "#include <m_append.h>" << '\n'
            << "#include <m_conversion.h>" << '\n'
            << "#include <m_copy.h>" << '\n'
            << "#include <m_gridded_fields.h>" << '\n'
            << "#include <m_xml.h>" << '\n'
            << "#include <m_select.h>" << '\n'
            << "#include <m_reduce.h>" << '\n'
            << "#include <m_nc.h>" << '\n'
            << "#include <m_delete.h>" << '\n'
            << "#include <m_extract.h>" << '\n'
            << "#include <m_ignore.h>" << '\n'
            << '\n'
            << '\n';
            
  std::cout << "namespace ARTS::Var {\n";
  for (auto& x : artsname.varname_group) {
    std::cout << "/*! " << x.second.varname_desc << '\n';
    std::cout << "@param[in,out] Workspace ws - An ARTS workspace\n@return "
                 "Reference to this workspace variable\n"
                 "*/\n";
    std::cout << "[[nodiscard]] inline ";
    std::cout << x.second.varname_group << '&' << ' ' << x.first
              << "(Workspace& ws) "
                 "noexcept { "
                 "return *static_cast<"
              << x.second.varname_group << " *>(ws[" << x.second.artspos
              << "]); "
                 "}\n\n";
  }
  std::cout << "}  // ARTS::Var \n\n";

  std::cout << "namespace ARTS::AgendaVar {\n";
  for (auto& x : artsname.group) {
    if (x.first == "Any") continue;
    
    std::cout << "class Workspace" << x.first << ' ' << '{' << '\n';
    std::cout << "  using type = " << x.first << ";\n";
    std::cout << "  std::size_t p;\n";
    std::cout << "  type* v;\n";
    std::cout << "public:\n";
    std::cout << "  Workspace" << x.first
              << "() noexcept : p(std::numeric_limits<std::size_t>::max()), "
                 "v(nullptr) {}\n";
    std::cout << "  Workspace" << x.first
              << "(std::size_t i, void * x) noexcept : p(i), "
                 "v(static_cast<type *>(x)) {}\n";
    std::cout << "  type& value() noexcept {return *v;}\n";
    std::cout << "  const type& value() const noexcept {return *v;}\n";
    std::cout
        << "  Workspace" << x.first
        << "& operator=(const type& t) noexcept {value() = t; return *this;}\n";
    std::cout << "  std::size_t pos() const noexcept {return p;}\n";
    std::cout << "  bool isnull() const noexcept {return v == nullptr;}\n";
    std::cout << '}' << ';' << '\n' << '\n';
  }
  for (auto& x : artsname.varname_group) {
    std::cout << "/*! " << x.second.varname_desc << '\n';
    std::cout << "@param[in,out] Workspace ws - An ARTS workspace\n";
    std::cout << "@return A class with a pointer to this variable and its "
                 "position in the workspace\n*/\n";
    std::cout << "[[nodiscard]] inline ";
    std::cout << "Workspace" << x.second.varname_group << ' ' << x.first
              << "(Workspace& ws) "
                 "noexcept { "
                 "return {"
              << x.second.artspos << ", ws[" << x.second.artspos
              << "]}; "
                 "}\n\n";
  }
  for (auto& x : artsname.group) {
    if (x.first == "Any") continue;
    
    std::cout << "/*! Creates in, and returns from, Workspace a/an " << x.first
              << '\n'
              << '\n';
    std::cout << "@param[in,out] Workspace ws - An ARTS workspace\n";
    std::cout << "@param[in] " << x.first
              << " inval - The default value the variable will have in "
                 "the workspace\n";
    std::cout << "@param[in] String name - The name the variable will have in "
                 "the workspace\n";
    std::cout << "@param[in] String desc - The description the variable will "
                 "have in the workspace (default: \"nodescription\")\n";
    std::cout << "@return A class with a pointer to this variable and its new "
                 "position in the workspace\n";
    std::cout << "*/\n";
    std::cout << "[[nodiscard]] inline\n";
    std::cout << "Workspace" << x.first << ' ' << x.first
              << "Create(\n            Workspace& ws,\n            const "
              << x.first
              << "& inval,\n            const String& name,\n            const "
                 "String& "
                 "desc=\"nodescription\") {\n";
    std::cout << "  const std::size_t ind = "
                 "std::size_t(ws.add_wsv_inplace({name.c_str(), desc.c_str(), "
              << x.second << "}));\n";
    std::cout << "  Workspace" << x.first << ' ' << "val{ind, ws[ind]};\n";
    std::cout << "  return val = inval;\n"
              << "}\n\n";
  }
  std::cout << "}  // ARTS::AgendaVar \n\n";

  std::cout << "namespace ARTS::Method {\n";

  for (auto& x : artsname.methodname_method) {
    // Skip methods using verbosity and Agenda methods (for now)
    if (x.agenda_method) continue;

    // Also skip create methods since these must be called via AgendaVar
    if (std::any_of(artsname.group.cbegin(), artsname.group.cend(),
                    [metname = x.name](auto& y) {
                      return (y.first + String("Create")) == metname;
                    }))
      continue;

    // Describe the method
    std::cout << "/*! " << x.desc << '\n';
    for (auto a : x.authors) std::cout << "@author " << a << '\n';
    std::cout << "\n"
                 "@param[in,out] Workspace ws - An ARTS workspace\n";
    for (std::size_t i = 0; i < x.gout.name.size(); i++)
      std::cout << "@param[out] " << x.gout.name[i] << " - " << x.gout.desc[i]
                << "\n";
    for (std::size_t i = 0; i < x.gin.name.size(); i++) {
      std::cout << "@param[in] " << x.gin.name[i] << " - " << x.gin.desc[i];
      if (x.gin.hasdefs[i]) std::cout << " (default: " << x.gin.defs[i] << ")";
      std::cout << '\n';
    }
    std::cout << "\nUse the ARTS documentation to read more on how the "
                 "workspace is manipulated\n";
    std::cout << "This interface function has been automatically generated\n";
    std::cout << "*/" << '\n';

    // Make the function
    std::cout << "inline void " << x.name << "(\n            Workspace& ws";

    // First put all GOUT variables
    if (x.pass_wsv_names) {
      for (std::size_t i = 0; i < x.gout.group.size(); i++) {
        std::cout << ',' << "\n            "
                  << "std::pair<" << x.gout.group[i] << ", String>" << '&'
                  << ' ' << x.gout.name[i];
      }
    } else {
      for (std::size_t i = 0; i < x.gout.group.size(); i++) {
        std::cout << ',' << "\n            " << x.gout.group[i] << '&' << ' '
                  << x.gout.name[i];
      }
    }

    // Second put all GIN variables that have no default argument
    if (x.pass_wsv_names) {
      for (std::size_t i = 0; i < x.gin.group.size(); i++) {
        if (not x.gin.hasdefs[i]) {
          std::cout << ',' << "\n            "
                    << "const std::pair<" << x.gin.group[i] << ", String>"
                    << '&' << ' ' << x.gin.name[i];
        }
      }
    } else {
      for (std::size_t i = 0; i < x.gin.group.size(); i++) {
        if (not x.gin.hasdefs[i]) {
          std::cout << ',' << "\n            "
                    << "const " << x.gin.group[i] << '&' << ' '
                    << x.gin.name[i];
        }
      }
    }

    // Lastly put all GIN variables that have a default argument
    if (x.pass_wsv_names) {
      for (std::size_t i = 0; i < x.gin.group.size(); i++) {
        if (x.gin.hasdefs[i]) {
          std::cout << ',' << "\n            "
                    << "const std::pair<" << x.gin.group[i] << ", String>"
                    << '&' << ' ' << x.gin.name[i] << '=' << '{'
                    << x.gin.defs[i] << ", \"" << x.gin.name[i] << "\"}";
        }
      }
    } else {
      for (std::size_t i = 0; i < x.gin.group.size(); i++) {
        if (x.gin.hasdefs[i]) {
          std::cout << ',' << "\n            "
                    << "const " << x.gin.group[i] << '&' << ' ' << x.gin.name[i]
                    << '=' << x.gin.defs[i];
        }
      }
    }

    // End of function definition and open function block
    std::cout << ')' << ' ' << '{' << '\n' << ' ' << ' ';

    // Call the ARTS auto_md.h function
    std::cout << x.name << '(';

    // We need the workspace if we input an Agenda or simply pass the workspace
    bool has_any = false;
    if (x.pass_workspace or x.agenda_method or
        std::any_of(
            x.gin.group.cbegin(), x.gin.group.cend(),
            [](auto& g) { return g == "Agenda" or g == "ArrayOfAgenda"; }) or
        std::any_of(x.in.varname.cbegin(), x.in.varname.cend(), [&](auto& g) {
          return artsname.varname_group.at(g).varname_group == "Agenda" or
                 artsname.varname_group.at(g).varname_group == "ArrayOfAgenda";
        })) {
      std::cout << "ws";
      has_any = true;
    }

    // First are all the outputs
    for (std::size_t i = 0; i < x.out.varname.size(); i++) {
      if (has_any) std::cout << ',' << ' ';
      has_any = true;
      std::cout << "Var::" << x.out.varname[i] << "(ws)";
    }

    // Second comes all the generic outputs
    for (std::size_t i = 0; i < x.gout.name.size(); i++) {
      if (has_any) std::cout << ',' << ' ';
      has_any = true;
      std::cout << x.gout.name[i];
      if (x.pass_wsv_names) std::cout << ".first";
    }

    // And their filenames if relevant
    if (x.pass_wsv_names) {
      for (std::size_t i = 0; i < x.gout.name.size(); i++) {
        if (has_any) std::cout << ',' << ' ';
        has_any = true;
        std::cout << x.gout.name[i] << ".second";
      }
    }

    // Then come all the inputs that are not also outputs
    for (std::size_t i = 0; i < x.in.varname.size(); i++) {
      if (std::any_of(
              x.out.varname.cbegin(), x.out.varname.cend(),
              [in = x.in.varname[i]](const auto& out) { return in == out; }))
        continue;
      if (has_any) std::cout << ',' << ' ';
      has_any = true;
      std::cout << "Var::" << x.in.varname[i] << "(ws)";
    }

    // Lastly are all the generic inputs, which cannot also be outputs
    for (std::size_t i = 0; i < x.gin.name.size(); i++) {
      if (has_any) std::cout << ',' << ' ';
      has_any = true;
      std::cout << x.gin.name[i];
      if (x.pass_wsv_names) std::cout << ".first";
    }

    // And their filenames if relevant
    if (x.pass_wsv_names) {
      for (std::size_t i = 0; i < x.gin.name.size(); i++) {
        if (has_any) std::cout << ',' << ' ';
        has_any = true;
        std::cout << x.gin.name[i] << ".second";
      }
    }

    // Check verbosity
    const bool has_verbosity =
        std::any_of(x.in.varname.cbegin(), x.in.varname.cend(),
                    [](auto& name) { return name == "verbosity"; });

    // Add verbosity of it does not exist
    if (not has_verbosity) {
      if (has_any) std::cout << ',' << ' ';
      has_any = true;
      std::cout << "Var::verbosity(ws)";
    }

    // Close the function call and the function itself
    std::cout << ')' << ';' << '\n' << '}' << '\n' << '\n' << '\n';
  }
  std::cout << "}  // ARTS::Method \n\n";

  std::cout << "namespace ARTS::AgendaMethod  {\n";

  for (auto& x : artsname.methodname_method) {
    // Skip methods using verbosity and Agenda methods (for now)
    if (x.agenda_method) continue;

    // Also skip create methods since these must be called via AgendaVar
    if (std::any_of(artsname.group.cbegin(), artsname.group.cend(),
                    [metname = x.name](auto& y) {
                      return (y.first + String("Create")) == metname;
                    }))
      continue;

    // Describe the method
    std::cout << "/*! " << x.desc << '\n';
    for (auto a : x.authors) std::cout << "@author " << a << '\n';
    std::cout << "\n"
                 "@param[in,out] Workspace ws - An ARTS workspace\n";
    for (std::size_t i = 0; i < x.gout.name.size(); i++)
      std::cout << "@param[out] " << x.gout.name[i] << " - " << x.gout.desc[i]
                << "\n";
    for (std::size_t i = 0; i < x.gin.name.size(); i++) {
      std::cout << "@param[in] " << x.gin.name[i] << " - " << x.gin.desc[i];
      if (x.gin.hasdefs[i]) std::cout << " (default: " << x.gin.defs[i] << ")";
      std::cout << '\n';
    }
    std::cout << "\nUse the ARTS documentation to read more on how the "
                 "workspace is manipulated\n";
    std::cout << "This interface function has been automatically generated\n";
    std::cout << "\n"
              << "@return MRecord to call this method\n";
    std::cout << "*/" << '\n';

    // Make the function
    std::cout << "[[nodiscard]] inline\nMRecord " << x.name
              << "(\n            [[maybe_unused]] Workspace& ws";

    // Check if we have the first input
    for (std::size_t i = 0; i < x.gout.group.size(); i++) {
      std::cout << ',' << '\n';
      std::cout << "                             AgendaVar::Workspace"
                << x.gout.group[i] << ' ' << x.gout.name[i];
    }

    // Second put all GIN variables that have no default argument
    for (std::size_t i = 0; i < x.gin.group.size(); i++) {
      if (not x.gin.hasdefs[i]) {
        std::cout << ',' << "\n";
        std::cout << "                       const AgendaVar::Workspace"
                  << x.gin.group[i] << ' ' << x.gin.name[i];
      }
    }

    // Lastly put all GIN variables that have a default argument
    for (std::size_t i = 0; i < x.gin.group.size(); i++) {
      if (x.gin.hasdefs[i]) {
        std::cout << ',' << "\n";
        std::cout << "                       const AgendaVar::Workspace"
                  << x.gin.group[i] << '&' << ' ' << x.gin.name[i] << '='
                  << "{}";
      }
    }

    // End of function definition and open function block
    std::cout << ')' << ' ' << '{' << '\n';

    for (std::size_t i = 0; i < x.gin.group.size(); i++) {
      if (x.gin.hasdefs[i]) {
        std::cout
            << "  static const auto " << x.gin.name[i]
            << "_default = AgendaVar::" << x.gin.group[i] << "Create(ws, "
            << x.gin.defs[i] << ",\n    \"" << x.name << '_' << x.gin.name[i]
            << "_autodefault"
            << "\",\n    \"auto generated variable with default from method "
               "definition\");\n";
      }
    }

    // Call the ARTS auto_md.h function
    std::cout << "  return MRecord(" << x.pos << ',' << ' '
              << "\n    ArrayOfIndex(" << '{';

    // First are all the outputs
    for (std::size_t i = 0; i < x.out.varpos.size(); i++) {
      std::cout << x.out.varpos[i] << ',' << ' ';
    }

    // Second comes all the generic outputs
    for (std::size_t i = 0; i < x.gout.name.size(); i++) {
      std::cout << "Index(" << x.gout.name[i] << ".pos())" << ',' << ' ';
    }
    std::cout << '}' << ')' << ',' << ' ' << "\n    ArrayOfIndex(" << '{';

    // Then come all the inputs that are not also outputs
    for (std::size_t i = 0; i < x.in.varpos.size(); i++) {
      std::cout << x.in.varpos[i] << ',' << ' ';
    }

    // Lastly are all the generic inputs, which cannot also be outputs
    for (std::size_t i = 0; i < x.gin.name.size(); i++) {
      if (x.gin.hasdefs[i])
        std::cout << x.gin.name[i] << ".isnull() ? Index(" << x.gin.name[i]
                  << "_default.pos()) : ";
      std::cout << "Index(" << x.gin.name[i] << ".pos())" << ',' << ' ';
    }

    std::cout << '}' << ')' << ',' << ' ';

    if (x.set_method)
      std::cout << "\n    TokVal{" << x.gin.name[0] << ".value()}";
    else
      std::cout << "\n    TokVal{}";

    std::cout << ", Agenda{}";

    // Close the function call and the function itself
    std::cout << ')' << ';' << '\n' << '}' << '\n' << '\n' << '\n';
  }
  std::cout << "}  // ARTS::AgendaMethod \n\n";

  std::cout << "namespace ARTS::AgendaExecute { \n\n";
  for (auto& x : artsname.agendaname_agenda) {
    std::cout << "/*! " << x.second.desc << '\n'
              << "@param[in,out] Workspace ws - An ARTS workspace\n"
              << "*/\n"
              << "inline void " << x.first << "(Workspace& ws) {\n  " << x.first
              << "Execute(ws";
    for (auto& name : x.second.outs) {
      std::cout << ',' << "\n            " << ' ' << "Var::" << name << "(ws)";
    }
    for (auto& name : x.second.ins) {
      if (not std::any_of(x.second.outs.cbegin(), x.second.outs.cend(),
                          [name](auto& outname) { return name == outname; }))
        std::cout << ',' << "\n            " << ' ' << "Var::" << name
                  << "(ws)";
    }
    std::cout << ",\n             Var::" << x.first << "(ws));\n}\n\n";
  }
  std::cout << "}  // ARTS::AgendaExecute \n\n";

  std::cout << "namespace ARTS::AgendaDefine { \n";
  std::cout << "/*! Append Records to an agenda */\n";
  std::cout << "template <typename ... Records>\nvoid Append(Agenda& ag, "
               "Records ... records) {\n"
            << "  for (auto& x: { MRecord(records)... })\n    "
               "ag.push_back(x);\n}\n\n";
  for (auto& x : artsname.agendaname_agenda) {
    if (artsname.varname_group.at(x.first).varname_group == "ArrayOfAgenda")
      continue;
    std::cout << "/*! " << x.second.desc << '\n'
              << "@param[in,out] Workspace ws - An ARTS workspace\n"
              << "@param[in] MRecords records - Any number of ARTS methods "
                 "from ARTS::AgendaMethod\n"
              << "*/\n"
              << "template <typename ... Records> "
              << "inline\nvoid " << x.first
              << "(Workspace& ws, Records ... records) {\n"
              << "  ARTS::Var::" << x.first << "(ws).resize(0);\n"
              << "  ARTS::Var::" << x.first << "(ws).set_name(\"" << x.first
              << "\");\n"
              << "  Append(ARTS::Var::" << x.first << "(ws), records...);"
              << "\n"
              << "  Var::" << x.first << "(ws).set_main_agenda();\n"
              << "  Var::" << x.first
              << "(ws).check(ws, Var::verbosity(ws));\n}\n\n";
  }
  std::cout << "}  // ARTS::AgendaDefine \n\n";
  
  // Make the main "startup"
  std::cout << "namespace ARTS {\n";
  std::cout <<
    "inline Workspace init(std::size_t screen=0, std::size_t file=0, std::size_t agenda=0) {\n"
    "  define_wsv_group_names();\n"
    "  Workspace::define_wsv_data();\n"
    "  Workspace::define_wsv_map();\n"
    "  define_md_data_raw();\n"
    "  expand_md_data_raw_to_md_data();\n"
    "  define_md_map();\n"
    "  define_agenda_data();\n"
    "  define_agenda_map();\n"
    "  define_species_data();\n"
    "  define_species_map();\n"
    "  global_data::workspace_memory_handler.initialize();\n"
    "\n"
    "  Workspace ws;\n"
    "  ws.initialize();\n"
    "  Var::verbosity(ws).set_screen_verbosity(screen);\n"
    "  Var::verbosity(ws).set_agenda_verbosity(agenda);\n"
    "  Var::verbosity(ws).set_file_verbosity(file);\n"
    "  Var::verbosity(ws).set_main_agenda(1);\n"
    "\n"
    "  #ifndef NDEBUG\n"
    "  ws.context = \"\";\n"
  "  #endif\n"
  "\n"
  "  return ws;"
  "}\n";
  std::cout << "}  // namespace::ARTS\n\n";

  std::cout << "#endif  // autoarts_h\n\n";
}
