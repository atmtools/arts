#include <auto_md.h>
#include <global_data.h>

#include <algorithm>
#include <filesystem>
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

void print_include_and_external(std::ofstream& core, std::ofstream& full) {
  core << "#include <auto_md.h>" << '\n' << "#include <global_data.h>" << '\n'
       << "#include <arts.h>" << '\n'
       << '\n'
       << '\n';
    
  full 
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
       << "#include <isotopologues.h>" << '\n'
       << '\n'
       << '\n';

  full << "extern String out_basename;\n\n";
}

void print_groups_and_namespaces(std::ofstream& core, const NameMaps& artsname) {
  core << "/*! An ARTS Workspace, the main class of ARTS */\n";
  core << "using Workspace=Workspace; \n\n";
  core << "/*! The ARTS constants namespace */\n";
  core << "namespace Constant { using namespace ::Constant; }\n\n";
  core << "/*! The ARTS conversions namespace */\n";
  core << "namespace Conversion { using namespace ::Conversion; }\n\n";
  core << "/*! The ARTS Group namespace */\n";
  core << "namespace Group {\n";
  for (auto& x : artsname.group) {
    if (x.first == "Any") continue;
    core << "/*! ARTS " << x.first << " type */\n";
    core << "using " << x.first << '=' << x.first << ';' << '\n' << '\n';
  }

  core << "/*! The ARTS Internal Groups namespace */\n";
  core << "namespace Internal {\n";
  core << "/*! A Tokenized Value.  Used purely in internal code */\n";
  core << "using TokVal=TokVal;\n" << '\n';
  core << "/*! A Method Record.  Used to pass methods to agendas */\n";
  core << "using MRecord=MRecord;\n" << '\n';
  core << "}  // namespace Internal\n";
  core << "}  // namespace Group \n\n";
}

void print_variables(std::ostream& core, const NameMaps& artsname) {
  core << "namespace Var {\n";
  for (auto& x : artsname.group) {
    if (x.first == "Any") continue;
   core << "/*! Workspace Variable class.  Used as default\n"
                 "  input to many Method and AgendaMethod.\n"
                 "  Note that it is not recommended to manually\n"
                 "  create this class as many methods makes distinct\n"
                 "  assumptions about the states that this class are\n"
                 "  allowed to be in */\n";
    core << "class " << x.first << ' ' << '{' << '\n';
    core << "  std::size_t p;\n";
    core << "  Group::" << x.first << "* v;\n";
    core << "public:\n";
    core << "  /*! Default construct.  DO NOT USE MANUALLY.  Leaves "
                 "islast() true and isnull() true */\n";
    core << "  " << x.first
              << "() noexcept : p(std::numeric_limits<std::size_t>::max()), "
                 "v(nullptr) {}\n\n";
    core << "  /*! Construct from existing Workspace.  DO NOT USE "
                 "MANUALLY. Leaves islast() false and isnull() false */\n";
    core << "  " << x.first
              << "(std::size_t i, void * x) noexcept : p(i), "
                 "v(static_cast<Group::" << x.first << " *>(x)) {}\n\n";
    core << "  /*! Construct from Group::" << x.first << ".  Leaves islast() true and "
                 "isnull() false */\n";
    core << "  " << x.first
              << "(const Group::" << x.first << "& val) noexcept : "
                 "p(std::numeric_limits<std::size_t>::max()), v(new Group::" << x.first << "(val)) "
                 "{}\n\n";
    core << "  /*! Delete data only when islast() is true and isnull() is "
                 "false */\n";
    core << "  ~" << x.first
              << "() noexcept {if (islast() and not isnull()) delete v;}\n\n";
    core
        << "  /*! Get value ref as Group::" << x.first << ".  Works when isnull() is false */\n";
    core << "  Group::" << x.first << "& value() noexcept {return *v;}\n\n";
    core
        << "  /*! Get value const-ref as Group::" << x.first << ".  Works when isnull() is false */\n";
    core << "  const Group::" << x.first << "& value() const noexcept {return *v;}\n\n";
    core << "  /*! Set value from Group::" << x.first << ".  Works when isnull() is "
                 "false */\n";
    core << "  " << x.first
              << "& operator=(const Group::" << x.first << "& t) noexcept {value() = t; return "
                 "*this;}\n\n";
    core
        << "  /*! Return the position of the variable in the Workspace */\n";
    core << "  std::size_t pos() const noexcept {return p;}\n\n";
    core << "  /*! Return true if there is no data */\n";
    core << "  bool isnull() const noexcept {return v == nullptr;}\n\n";
    core << "  /*! Return true if data is not in the Workspace */\n";
    core << "  bool islast() const noexcept {return p == "
                 "std::numeric_limits<std::size_t>::max();}\n\n";
    core << "  /*! Name of variable.  Must be in the workspace */\n";
    core << "  const Group::String& name() const noexcept {return "
                 "Workspace::wsv_data[p].Name();}\n\n";

    // NOTE:  Don't add more groups here.  The ones that are here are ugly
    // enough as it is.  Just define an output operator...
    if (x.first not_eq "Ppath" and x.first not_eq "TessemNN" and
        x.first not_eq "Timer" and x.first not_eq "ArrayOfPpath") {
      core << "  /*! Output to stream of internal variable */\n";
      core << "  friend std::ostream& operator<<(std::ostream& os, const "
                << x.first
                << "& var) {if (var.isnull()) return os << \"NULLDATA\"; else "
                   "return os << var.value();}\n";
    }

    core << '}' << ';' << '\n' << '\n';
  }
  for (auto& x : artsname.varname_group) {
    core << "/*! " << x.second.varname_desc << '\n';
    core << "@param[in,out] Workspace ws - An ARTS workspace\n";
    core << "@return A class with a pointer to this variable and its "
                 "position in the workspace\n*/\n";
    core << "[[nodiscard]] inline ";
    core << x.second.varname_group << ' ' << x.first
              << "(Workspace& ws) "
                 "noexcept { "
                 "return {"
              << x.second.artspos << ", ws[" << x.second.artspos
              << "]}; "
                 "}\n\n";
  }
  for (auto& x : artsname.group) {
    if (x.first == "Any") continue;

    core << "/*! Creates in, and returns from, Workspace a/an " << x.first
              << '\n'
              << '\n';
    core << "@param[in,out] Workspace ws - An ARTS workspace\n";
    core << "@param[in] " << x.first
              << " inval - The default value the variable will have in "
                 "the workspace\n";
    core << "@param[in] String name - The name the variable will have in "
                 "the workspace\n";
    core << "@param[in] String desc - The description the variable will "
                 "have in the workspace (default: \"nodescription\")\n";
    core << "@return A class with a pointer to this variable and its new "
                 "position in the workspace\n";
    core << "*/\n";
    core << "[[nodiscard]] inline\n";
    core
        << x.first << ' ' << x.first
        << "Create(\n            Workspace& ws,\n            const Group::"
        << x.first
        << "& inval,\n            const Group::String& name,\n            "
           "const Group::"
           "String& "
           "desc=\"nodescription\") {\n";
    core << "  const std::size_t ind = "
                 "std::size_t(ws.add_wsv_inplace({name.c_str(), desc.c_str(), "
              << x.second << "}));\n";
    core << "  " << x.first << ' ' << "val{ind, ws[ind]};\n";
    core << "  return val = inval;\n"
              << "}\n\n";
  }
  core << "}  // namespace Var \n\n";
}

std::string spaces(std::size_t n) {
  std::string s = "";
  for (std::size_t i=0; i<n; i++)
    s += " ";
  return s;
}

void print_gin_methods(std::ofstream& full, const NameMaps& artsname) {
  for (auto& x : artsname.methodname_method) {
    // Skip methods using verbosity and Agenda methods (for now)
    if (x.agenda_method) continue;

    // Also skip create methods since these must be called via Var
    if (std::any_of(artsname.group.cbegin(), artsname.group.cend(),
                    [metname = x.name](auto& y) {
                      return (y.first + String("Create")) == metname;
                    }))
      continue;

    // Describe the method
    full << "/*! " << x.desc << '\n';
    for (auto a : x.authors) full << "@author " << a << '\n';
    full << "\n"
                 "@param[in,out] Workspace ws - An ARTS workspace\n";
    for (std::size_t i = 0; i < x.gout.name.size(); i++)
      full << "@param[out] " << x.gout.name[i] << " - " << x.gout.desc[i]
                << "\n";
    for (std::size_t i = 0; i < x.gin.name.size(); i++) {
      full << "@param[in] " << x.gin.name[i] << " - " << x.gin.desc[i];
      if (x.gin.hasdefs[i]) full << " (default: " << x.gin.defs[i] << ")";
      full << '\n';
    }
    full << "\nUse the ARTS documentation to read more on how the "
                 "workspace is manipulated\n";
    full << "This interface function has been automatically generated\n";
    full << "*/" << '\n';

    // Make the function
    full << "inline void " << x.name << "(Workspace& ws";

    // First put all GOUT variables
    for (std::size_t i = 0; i < x.gout.group.size(); i++) {
      full << ',' << '\n' << spaces(x.name.size() + 13) << "Var::" << x.gout.group[i] << ' '
                << x.gout.name[i];
    }

    // Second put all GIN variables that have no default argument
    for (std::size_t i = 0; i < x.gin.group.size(); i++) {
      if (not x.gin.hasdefs[i]) {
        full << ',' << '\n' << spaces(x.name.size() + 13 - 6)
                  << "const Var::" << x.gin.group[i] << ' ' << x.gin.name[i];
      }
    }

    // Lastly put all GIN variables that have a default argument
    for (std::size_t i = 0; i < x.gin.group.size(); i++) {
      if (x.gin.hasdefs[i]) {
        if (x.gin.defs[i] == "{}") {
          full << ',' << '\n' << spaces(x.name.size() + 13 - 6)
                    << "const Var::" << x.gin.group[i] << ' ' << x.gin.name[i]
                    << '=' << "Group::" << x.gin.group[i] << x.gin.defs[i];
        } else {
          full << ',' << '\n' << spaces(x.name.size() + 13 - 6)
                    << "const Var::" << x.gin.group[i] << ' ' << x.gin.name[i]
                    << '=' << "Group::" << x.gin.group[i] << '{'
                    << x.gin.defs[i] << '}';
        }
      }
    }

    // End of function definition and open function block
    full << ')' << ' ' << '{' << '\n';

    // Output variables have to be on the Workspace
    if (x.gout.group.size()) full << ' ';
    for (std::size_t i = 0; i < x.gout.group.size(); i++) {
      full << " if (" << x.gout.name[i]
                << ".islast()) {\n    throw std::runtime_error(\""
                << x.gout.name[i] << " needs to be a defined "
                << x.gout.group[i] << " since it is output of " << x.name
                << "\");\n  }";
    }
    if (x.gout.group.size()) full << '\n' << '\n';

    // Call the ARTS auto_md.h function
    full << ' ' << ' ' << x.name << '(';

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
      full << "ws";
      has_any = true;
    }

    // First are all the outputs
    for (std::size_t i = 0; i < x.out.varname.size(); i++) {
      if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
      has_any = true;
      full << "Var::" << x.out.varname[i] << "(ws).value()";
    }

    // Second comes all the generic outputs
    for (std::size_t i = 0; i < x.gout.name.size(); i++) {
      if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
      has_any = true;
      full << x.gout.name[i];
      full << ".value()";
    }

    // And their filenames if relevant
    if (x.pass_wsv_names) {
      for (std::size_t i = 0; i < x.gout.name.size(); i++) {
        if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
        has_any = true;
        full << x.gout.name[i] << ".name()";
      }
    }

    // Then come all the inputs that are not also outputs
    for (std::size_t i = 0; i < x.in.varname.size(); i++) {
      if (std::any_of(
              x.out.varname.cbegin(), x.out.varname.cend(),
              [in = x.in.varname[i]](const auto& out) { return in == out; }))
        continue;
      if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
      has_any = true;
      full << "Var::" << x.in.varname[i] << "(ws).value()";
    }

    // Lastly are all the generic inputs, which cannot also be outputs
    for (std::size_t i = 0; i < x.gin.name.size(); i++) {
      if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
      has_any = true;
      full << x.gin.name[i];
      full << ".value()";
    }

    // And their filenames if relevant
    if (x.pass_wsv_names) {
      for (std::size_t i = 0; i < x.gin.name.size(); i++) {
        if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
        has_any = true;
        full << x.gin.name[i] << ".islast() ? Group::String{\""
                  << x.gin.name[i] << "\"} : " << x.gin.name[i] << ".name()";
      }
    }

    // Check verbosity
    const bool has_verbosity =
        std::any_of(x.in.varname.cbegin(), x.in.varname.cend(),
                    [](auto& name) { return name == "verbosity"; });

    // Add verbosity of it does not exist
    if (not has_verbosity) {
      if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
      has_any = true;
      full << "Var::verbosity(ws).value()";
    }

    // Close the function call and the function itself
    full << ')' << ';' << '\n' << '}' << '\n' << '\n' << '\n';
  }
}

void print_full_methods(std::ofstream& full, const NameMaps& artsname) {
  for (auto& x : artsname.methodname_method) {
    // Skip methods using verbosity and Agenda methods (for now)
    if (x.agenda_method) continue;

    // Also skip create methods since these must be called via Var
    if (std::any_of(artsname.group.cbegin(), artsname.group.cend(),
                    [metname = x.name](auto& y) {
                      return (y.first + String("Create")) == metname;
                    }))
      continue;

    // Skip the 'silly' methods
    if (not x.out.varname.size() and not x.gout.name.size() and
        not x.in.varname.size() and not x.gin.name.size())
      continue;

    // Describe the method
    full << "/*! " << x.desc << '\n';
    for (auto a : x.authors) full << "@author " << a << '\n';
    full << "\n"
                 "@param[in,out] Workspace ws - An ARTS workspace\n";
    for (std::size_t i = 0; i < x.out.varname.size(); i++) {
      if (std::any_of(
              x.in.varname.cbegin(), x.in.varname.cend(),
              [out = x.out.varname[i]](const auto& in) { return in == out; }))
        full << "@param[in,out] ";
      else
        full << "@param[out] ";
      full << x.out.varname[i] << " - as *Var::" << x.out.varname[i]
                << "(ws)*\n";
    }
    for (std::size_t i = 0; i < x.gout.name.size(); i++)
      full << "@param[out] " << x.gout.name[i] << " - " << x.gout.desc[i]
                << "\n";
    for (std::size_t i = 0; i < x.in.varname.size(); i++) {
      if (std::any_of(
              x.out.varname.cbegin(), x.out.varname.cend(),
              [in = x.in.varname[i]](const auto& out) { return in == out; }))
        continue;
      full << "@param[in] " << x.in.varname[i]
                << " - as *Var::" << x.in.varname[i] << "(ws)*\n";
    }
    for (std::size_t i = 0; i < x.gin.name.size(); i++) {
      full << "@param[in] " << x.gin.name[i] << " - " << x.gin.desc[i];
      if (x.gin.hasdefs[i]) full << " (default: " << x.gin.defs[i] << ")";
      full << '\n';
    }
    full << "\nUse the ARTS documentation to read more on how the "
                 "workspace is manipulated\n";
    full << "This interface function has been automatically generated\n";
    full << "*/" << '\n';

    // Make the function
    full << "inline void " << x.name
              << "(Workspace& ws [[maybe_unused]]";
    
    // First put all OUT variables
    for (std::size_t i = 0; i < x.out.varname.size(); i++)
      full << ',' << '\n' << spaces(x.name.size() + 13) << "Group::"
                << artsname.varname_group.at(x.out.varname[i]).varname_group
                << '&' << ' ' << x.out.varname[i];

    // Second put all GOUT variables
    for (std::size_t i = 0; i < x.gout.group.size(); i++)
      full << ',' << '\n' << spaces(x.name.size() + 13) << "Group::" << x.gout.group[i]
                << '&' << ' ' << x.gout.name[i];

    // Third put all the inputs that are not also outputs
    for (std::size_t i = 0; i < x.in.varname.size(); i++) {
      if (std::any_of(
              x.out.varname.cbegin(), x.out.varname.cend(),
              [in = x.in.varname[i]](const auto& out) { return in == out; }))
        continue;
      full << ',' << '\n' << spaces(x.name.size() + 13 - 6) << "const Group::"
                << artsname.varname_group.at(x.in.varname[i]).varname_group
                << '&' << ' ' << x.in.varname[i];
    }

    // Last put all GIN variables that have no default argument
    for (std::size_t i = 0; i < x.gin.group.size(); i++)
      full << ',' << '\n' << spaces(x.name.size() + 13 - 6) << "const Group::" << x.gin.group[i] << '&'
                << ' ' << x.gin.name[i];

    // End of function definition and open function block
    full << ')' << ' ' << '{' << '\n';

    // Call the ARTS auto_md.h function
    full << ' ' << ' ' << x.name << '(';

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
      full << "ws";
      has_any = true;
    }

    // First are all the outputs
    for (std::size_t i = 0; i < x.out.varname.size(); i++) {
      if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
      has_any = true;
      full << x.out.varname[i];
    }

    // Second comes all the generic outputs
    for (std::size_t i = 0; i < x.gout.name.size(); i++) {
      if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
      has_any = true;
      full << x.gout.name[i];
    }

    // And their filenames if relevant
    if (x.pass_wsv_names) {
      for (std::size_t i = 0; i < x.gout.name.size(); i++) {
        if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
        has_any = true;
        full << "Group::String{\"" << x.gout.name[i] << "\"}";
      }
    }

    // Then come all the inputs that are not also outputs
    for (std::size_t i = 0; i < x.in.varname.size(); i++) {
      if (std::any_of(
              x.out.varname.cbegin(), x.out.varname.cend(),
              [in = x.in.varname[i]](const auto& out) { return in == out; }))
        continue;
      if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
      has_any = true;
      full << x.in.varname[i];
    }

    // Lastly are all the generic inputs, which cannot also be outputs
    for (std::size_t i = 0; i < x.gin.name.size(); i++) {
      if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
      has_any = true;
      full << x.gin.name[i];
    }

    // And their filenames if relevant
    if (x.pass_wsv_names) {
      for (std::size_t i = 0; i < x.gin.name.size(); i++) {
        if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
        has_any = true;
        full << "Group::String{\"" << x.gin.name[i] << "\"}";
      }
    }

    // Check verbosity
    const bool has_verbosity =
        std::any_of(x.in.varname.cbegin(), x.in.varname.cend(),
                    [](auto& name) { return name == "verbosity"; });

    // Add verbosity of it does not exist
    if (not has_verbosity) {
      if (has_any) full << ',' << '\n' << spaces(x.name.size() + 2 + 1);
      has_any = true;
      full << "Var::verbosity(ws).value()";
    }

    // Close the function call and the function itself
    full << ')' << ';' << '\n' << '}' << '\n' << '\n' << '\n';
  }
}

void print_methods(std::ofstream& full, const NameMaps& artsname) {
  full << "#pragma clang diagnostic push\n";
  full << "#pragma clang diagnostic ignored \"-Wbraced-scalar-init\"\n";
  full << "namespace Method {\n";
  print_gin_methods(full, artsname);
  print_full_methods(full, artsname);
  full << "}  // namespace Method \n\n";
  full << "#pragma clang diagnostic pop\n";
}

void print_agenda_methods(std::ofstream& core, const NameMaps& artsname) {
  core << "namespace Method {\n";

  for (auto& x : artsname.methodname_method) {
    // Skip methods using verbosity and Agenda methods (for now)
    if (x.agenda_method) continue;

    // Also skip create methods since these must be called via Var
    if (std::any_of(artsname.group.cbegin(), artsname.group.cend(),
                    [metname = x.name](auto& y) {
                      return (y.first + String("Create")) == metname;
                    }))
      continue;

    // Describe the method
    core << "/*! " << x.desc << '\n';
    for (auto a : x.authors) core << "@author " << a << '\n';
    core << "\n"
                 "@param[in,out] Workspace ws - An ARTS workspace\n";
    for (std::size_t i = 0; i < x.gout.name.size(); i++)
      core << "@param[out] " << x.gout.name[i] << " - " << x.gout.desc[i]
                << "\n";
    for (std::size_t i = 0; i < x.gin.name.size(); i++) {
      core << "@param[in] " << x.gin.name[i] << " - " << x.gin.desc[i];
      if (x.gin.hasdefs[i]) core << " (default: " << x.gin.defs[i] << ")";
      core << '\n';
    }
    core << "\nUse the ARTS documentation to read more on how the "
                 "workspace is manipulated\n";
    core << "This interface function has been automatically generated\n";
    core << "\n"
              << "@return MRecord to call this method\n";
    core << "*/" << '\n';

    // Make the function
    core << "[[nodiscard]] inline\nGroup::Internal::MRecord " << x.name
              << "(Workspace& ws [[maybe_unused]]";

    // Check if we have the first input
    for (std::size_t i = 0; i < x.gout.group.size(); i++) {
      core << ',' << '\n' << spaces(x.name.size() + 26);
      core << "Var::" << x.gout.group[i] << ' '
                << x.gout.name[i];
    }

    // Second put all GIN variables that have no default argument
    for (std::size_t i = 0; i < x.gin.group.size(); i++) {
      if (not x.gin.hasdefs[i]) {
        core << ',' << '\n' << spaces(x.name.size() + 26 - 6);
        core << "const Var::" << x.gin.group[i] << ' '
                  << x.gin.name[i];
      }
    }

    // Lastly put all GIN variables that have a default argument
    for (std::size_t i = 0; i < x.gin.group.size(); i++) {
      if (x.gin.hasdefs[i]) {
        core << ',' << '\n' << spaces(x.name.size() + 26 - 6);
        core << "const Var::" << x.gin.group[i] << '&'
                  << ' ' << x.gin.name[i] << '=' << "{}";
      }
    }

    // End of function definition and open function block
    core << ')' << ' ' << '{' << '\n';

    // Output variables have to be on the Workspace
    if (x.gout.group.size() or x.gin.group.size()) core << ' ';
    for (std::size_t i = 0; i < x.gout.group.size(); i++) {
      core << " if (" << x.gout.name[i]
                << ".islast()) {\n    throw std::runtime_error(\""
                << x.gout.name[i] << " needs to be a defined Workspace"
                << x.gout.group[i] << " since it is output of " << x.name
                << "\");\n  }";
    }
    for (std::size_t i = 0; i < x.gin.group.size(); i++) {
      if (x.gin.hasdefs[i])
        core << " if (not " << x.gin.name[i] << ".isnull() and "
                  << x.gin.name[i]
                  << ".islast()) {\n    throw std::runtime_error(\""
                  << x.gin.name[i] << " needs to be a defined Workspace"
                  << x.gin.group[i]
                  << " (or left default) since it is agenda input to " << x.name
                  << "\");\n  }";
      else
        core << " if (" << x.gin.name[i]
                  << ".islast()) {\n    throw std::runtime_error(\""
                  << x.gin.name[i] << " needs to be a defined Workspace"
                  << x.gin.group[i] << " since it is agenda input to " << x.name
                  << "\");\n  }";
    }
    if (x.gout.group.size() or x.gin.group.size()) core << '\n' << '\n';

    for (std::size_t i = 0; i < x.gin.group.size(); i++) {
      if (x.gin.hasdefs[i]) {
        core
            << "  static const auto " << x.gin.name[i]
            << "_default = Var::" << x.gin.group[i] << "Create(ws, "
            << x.gin.defs[i] << ",\n    \"" << x.name << '_' << x.gin.name[i]
            << "_autodefault"
            << "\",\n    \"auto generated variable with default from method "
               "definition\");\n";
      }
    }

    // Call the ARTS auto_md.h function
    core << "  return Group::Internal::MRecord(" << x.pos << ',' << '\n'
              << spaces(34) << "Group::ArrayOfIndex(" << '{';

    // First are all the outputs
    for (std::size_t i = 0; i < x.out.varpos.size(); i++) {
      core << x.out.varpos[i] << ',' << '\n' << spaces(34 + 21);
    }

    // Second comes all the generic outputs
    for (std::size_t i = 0; i < x.gout.name.size(); i++) {
      core << "Group::Index(" << x.gout.name[i] << ".pos())" << ',' << '\n' << spaces(34 + 21);
    }
    core << '}' << ')' << ',' << '\n' << spaces(34) << "Group::ArrayOfIndex("
              << '{';

    // Then come all the inputs that are not also outputs
    for (std::size_t i = 0; i < x.in.varpos.size(); i++) {
      core << x.in.varpos[i] << ',' << '\n' << spaces(34 + 21);
    }

    // Lastly are all the generic inputs, which cannot also be outputs
    for (std::size_t i = 0; i < x.gin.name.size(); i++) {
      if (x.gin.hasdefs[i])
        core << x.gin.name[i] << ".isnull() ? Group::Index(" << x.gin.name[i]
                  << "_default.pos()) :" << '\n' << spaces(34 + 21 + 12 + x.gin.name[i].size());
      core << "Group::Index(" << x.gin.name[i] << ".pos())" << ',' << '\n' << spaces(34 + 21);
    }

    core << '}' << ')' << ',' << ' ';

    if (x.set_method)
      core << '\n' << spaces(34) << "Group::Internal::TokVal{" << x.gin.name[0]
                << ".value()}";
    else
      core << '\n' << spaces(34) << "Group::Internal::TokVal{}";

    core << ',' << '\n' << spaces(34) << "Group::Agenda{}";

    // Close the function call and the function itself
    core << ')' << ';' << '\n' << '}' << '\n' << '\n' << '\n';
  }
  core << "}  // namespace Method \n\n";
}

void print_agenda_execute(std::ofstream& core, const NameMaps& artsname) {
  core << "namespace Execute { \n\n";
  for (auto& x : artsname.agendaname_agenda) {
    core << "/*! " << x.second.desc << '\n'
              << "@param[in,out] Workspace ws - An ARTS workspace\n"
              << "*/\n"
              << "inline void " << x.first << "(Workspace& ws) {\n  " << x.first
              << "Execute(ws";
    for (auto& name : x.second.outs) {
      core << ',' << '\n' << spaces(x.first.size() + 2 + 1 + 7) << "Var::" << name
                << "(ws).value()";
    }
    for (auto& name : x.second.ins) {
      if (not std::any_of(x.second.outs.cbegin(), x.second.outs.cend(),
                          [name](auto& outname) { return name == outname; }))
        core << ',' << '\n' << spaces(x.first.size() + 2 + 1 + 7) << "Var::" << name
                  << "(ws).value()";
    }
    core << ',' << '\n' << spaces(x.first.size() + 2 + 1 + 7) << "Var::" << x.first << "(ws).value());\n}\n\n";
  }
  core << "}  // namespace Execute \n\n";
}

void print_agenda_define(std::ofstream& core, const NameMaps& artsname) {
  core << "namespace Define { \n";
  core << "/*! Append Records to an agenda */\n";
  core
      << "template <typename ... Records>\nvoid Append(Group::Agenda& ag, "
         "Records ... records) {\n"
      << "  for (auto& x: { Group::Internal::MRecord(records)... })\n    "
         "ag.push_back(x);\n}\n\n";
  for (auto& x : artsname.agendaname_agenda) {
    if (artsname.varname_group.at(x.first).varname_group == "ArrayOfAgenda")
      continue;
    core << "/*! " << x.second.desc << '\n'
              << "@param[in,out] Workspace ws - An ARTS workspace\n"
              << "@param[in] MRecords records - Any number of ARTS methods "
                 "from AgendaMethod\n"
              << "*/\n"
              << "template <typename ... Records> "
              << "\nvoid " << x.first
              << "(Workspace& ws, Records ... records) {\n"
              << "  Var::" << x.first << "(ws).value().resize(0);\n"
              << "  Var::" << x.first << "(ws).value().set_name(\"" << x.first
              << "\");\n"
              << "  Append(Var::" << x.first << "(ws).value(), records...);"
              << "\n"
              << "  Var::" << x.first
              << "(ws).value().check(ws, Var::verbosity(ws).value());\n}\n\n";
  }
  core << "}  // namespace Define \n\n";
}

void print_agendas(std::ofstream& core, const NameMaps& artsname) {
  core << "/*! ARTS Agenda interaction namespace\n\n   Will only be "
               "populated by namespaces.\n   Users of the API can define their "
               "Agendas\n   in the Agenda namespace\n*/\n";
  core << "namespace Agenda {\n";
  print_agenda_methods(core, artsname);
  print_agenda_execute(core, artsname);
  print_agenda_define(core, artsname);
  core << "}  // namespace Agenda\n";
}

void print_species_identification(std::ofstream& core) {
  core << "namespace Species {\n";
  core << "  using namespace ::Species;\n";
  core << "}  // namespace Species \n\n";
}

void print_startup(std::ofstream& full) {
  full << "/*! Create a Workspace and set its main verbosity\n\n"
            << "  @param[in] screen Screen verbosity\n"
            << "  @param[in] file File verbosity\n"
            << "  @param[in] agenda Agenda verbosity\n"
            << "  @param[in] basename Default basename for output variables\n"
            << "  @param[in] numthreads OpenMP thread count (defaults to max "
               "if invalid count)\n"
            << "  @return Workspace a full ARTS Workspace\n"
            << "*/\n";
            full
      << "inline Workspace init(std::size_t screen=0, std::size_t file=0, "
         "std::size_t agenda=0, const Group::String& basename=\"arts\", "
         "[[maybe_unused]] int numthreads=0) {\n"
#ifdef _OPENMP
         "  omp_set_num_threads(numthreads < 1 ? arts_omp_get_max_threads() : "
         "numthreads > arts_omp_get_max_threads() ? arts_omp_get_max_threads() "
         ": numthreads);\n"
         "\n"
#endif
         "  define_wsv_group_names();\n"
         "  Workspace::define_wsv_data();\n"
         "  Workspace::define_wsv_map();\n"
         "  define_md_data_raw();\n"
         "  expand_md_data_raw_to_md_data();\n"
         "  define_md_map();\n"
         "  define_agenda_data();\n"
         "  define_agenda_map();\n"
         "  global_data::workspace_memory_handler.initialize();\n"
         "\n"
         "  Workspace ws;\n"
         "  ws.initialize();\n"
         "  Var::verbosity(ws).value().set_screen_verbosity(screen);\n"
         "  Var::verbosity(ws).value().set_agenda_verbosity(agenda);\n"
         "  Var::verbosity(ws).value().set_file_verbosity(file);\n"
         "  Var::verbosity(ws).value().set_main_agenda(1);\n"
         "\n"
         "  out_basename = basename;\n"
         "\n"
         "#ifndef NDEBUG\n"
         "  ws.context = \"\";\n"
         "#endif\n"
         "\n"
         "  return ws;\n"
         "}\n";
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

  const std::filesystem::path dir(".");

  // Open file handlers
  std::ofstream core(dir / "autoarts_core.h");
  std::ofstream full(dir / "autoarts.h");

  const auto artsname = NameMaps();

  core << "#ifndef autoarts_core_h\n"
       << "#define autoarts_core_h\n"
       << '\n';

  full << "#ifndef autoarts_h\n"
       << "#define autoarts_h\n"
       << '\n'
       << "#include <autoarts_core.h>\n"
       << '\n';

  print_include_and_external(core, full);

  core << "namespace ARTS {\n";
  full << "namespace ARTS {\n";

  print_groups_and_namespaces(core, artsname);

  print_variables(core, artsname);

  print_methods(full, artsname);

  print_agendas(core, artsname);

  print_species_identification(core);

  print_startup(full);

  core << "}  // namespace::ARTS\n\n";
  full << "}  // namespace::ARTS\n\n";

  core << "#endif  // autoarts_core_h\n";
  full << "#endif  // autoarts_h\n";
}
