#include <auto_md.h>
#include <global_data.h>

#include <fstream>

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

void includes(std::ofstream& os) {
  os << "#include <auto_md.h>" << '\n';
  os << "#include <xml_io.h>" << '\n';
  os << "#include <py_macros.h>" << '\n';
  os << "#include <python_interface.h>" << '\n';
}

void workspace_variables(std::ofstream& os, const NameMaps& arts) {
  for (auto& [name, data]: arts.varname_group) {
    os << "  ws.def_property(\"" << name << "\",\n";
    os << "    [](Workspace& w) {return * reinterpret_cast<" << data.varname_group << " *> ws["<<data.artspos<<"];},\n";
    os << "    [](Workspace& w, " << data.varname_group << " val) {(* reinterpret_cast<" << data.varname_group << " *> ws["<<data.artspos<<"]) = std::move(val);},\n";
    os << "    py::doc(LR\"-VARNAME_DESC-("<<data.varname_desc<<")-VARNAME_DESC-\");";
    os << '\n' << '\n';
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

  std::ofstream py_workspace("py_workspace.cc");

  includes(py_workspace);

  const auto artsname = NameMaps();

  py_workspace << '\n';
  py_workspace << "namespace Python {\n";
  py_workspace << "void py_auto_workspace(py::class_<Workspace>& ws) {\n";
  workspace_variables(py_workspace, artsname);
  py_workspace << "}\n";
  py_workspace << "}  // namespace Python\n";
}
