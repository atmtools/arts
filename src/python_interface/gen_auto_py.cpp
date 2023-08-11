#include <workspace.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <ostream>
#include <sstream>

#include "pydocs.h"

std::string fix_type(const std::string& name) {
  if ("Numeric" == name or "Index" == name) return name + "_";
  return "std::shared_ptr<" + name + ">";
}

std::string variable(const std::string& name,
                     const WorkspaceVariableRecord& wsv) try {
  std::ostringstream os;

  os << "  ws.def_property(\"" << name << "\",";
  os << R"--(
    py::cpp_function([](Workspace& w) -> )--"
     << fix_type(wsv.type) << R"--( {
      return w.share(")--" << name << R"--(") -> share_unsafe<)--" << wsv.type<< R"--(>();
    }, py::return_value_policy::reference_internal, py::keep_alive<0, 1>()),
    [](Workspace& w, )--"
     << fix_type(wsv.type) << R"--( val) -> void {
      w.set(")--" << name << R"--(", std::make_shared<Wsv>(std::move(val)--";
      if (wsv.type == "Numeric" or wsv.type == "Index") os << ".val"; 
      os << R"--()));
)--";

  if (wsv.type == "Agenda")
    os << "      auto& ws_val = w.get<"<<wsv.type<<">(\""<<name<<"\");\n      ws_val.set_name(\"" << name
       << "\");\n      ws_val.finalize();\n";
  if (wsv.type == "ArrayOfAgenda")
    os << "      auto& ws_val = w.get<"<<wsv.type<<">(\""<<name<<"\");\n      for (auto&ws_value: ws_val) {\n        ws_value.set_name(\""
       << name << "\");\n        ws_value.finalize();\n      }\n";

  os << "    }, py::doc(R\"-x-(:class:`~pyarts.arts." << wsv.type << "` "
     << unwrap_stars(wsv.desc) << "\n\n";

  if (wsv.type == "Agenda" or wsv.type == "ArrayOfAgenda") {
    os << get_agenda_io(name);
  }

  if (wsv.default_value) {
    os << "\nDefault value\n"
          "-------------\n\n``"
       << to_defval_str(*wsv.default_value, wsv.type) << "``\n\n";
  }

  os << variable_used_by(name) << '\n';

  os << "\n\nGeneric workspace methods that can generate or use " << name
     << "\n"
     << String(51 + name.size(), '-') << "\n\nSee :class:`~pyarts.arts."
     << wsv.type << "`";
  if (wsv.type not_eq "Any") os << " and/or :class:`~pyarts.arts.Any`\n";

  return fix_newlines(os.str()) + ")-x-\"));\n\n";
} catch (std::exception& e) {
  std::cerr << "Error in variable " << std::quoted(name) << ":\n"
            << e.what() << '\n';
  std::exit(1);
}

void variables(const std::string& fname) {
  const auto& wsvs = workspace_variables();

  std::ofstream os(fname);

  os << R"--(#include <python_interface.h>

namespace Python {
void py_auto_wsv(artsclass<Workspace>& ws [[maybe_unused]]) {
)--";

  for (auto& [name, wsv] : wsvs) {
    os << variable(name, wsv) << std::flush;
  }
  os << "}\n}  // namespace Python\n";
}

std::vector<std::string> supergeneric_type(const std::string& t) {
  std::vector<std::string> types;
  std::string type = t;
  auto ptr = type.find(',');
  while (ptr != type.npos) {
    types.push_back(type.substr(0, ptr));
    type = type.substr(ptr + 1);
    ptr = type.find(',');
  }
  types.push_back(type);

  for (auto& T : types) {
    while (T.back() == ' ') T.pop_back();
    while (T.front() == ' ') T.erase(0, 1);
  }

  for (auto& T : types) {
    T = fix_type(T);
  }

  return types;
}

std::vector<std::string> unique_sorted(const std::vector<std::string>& v) {
  std::vector<std::string> out = v;
  std::ranges::sort(out);
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

std::vector<std::vector<std::string>> supergeneric_types(
    const WorkspaceMethodInternalRecord& wsm) {
  std::vector<std::vector<std::string>> out;

  for (auto& t : wsm.gout_type) {
    if (t.find(',') == t.npos) continue;

    const auto types = supergeneric_type(t);

    if (not out.size()) {
      out.resize(types.size());
    }
    for (std::size_t i = 0; i < types.size(); i++) {
      out[i].push_back(types[i]);
    }
  }

  for (auto& t : wsm.gin_type) {
    if (t.find(',') == t.npos) continue;

    const auto types = supergeneric_type(t);

    if (not out.size()) {
      out.resize(types.size());
    }
    for (std::size_t i = 0; i < types.size(); i++) {
      out[i].push_back(types[i]);
    }
  }
  return out;
}

std::string method_g_types(const std::string& type) {
  if (type == "Any") {
    return "PyWsvValue";
  }

  if (type.find(',') != type.npos) {
    std::ostringstream os;
    os << "std::variant<";
    bool first = true;
    for (auto&& t : unique_sorted(supergeneric_type(type))) {
      if (not first) os << ", ";
      first = false;
      os << t;
    }
    os << ">";
    return os.str();
  }

  return fix_type(type);
}

std::string method_arguments(const WorkspaceMethodInternalRecord& wsm) {
  static const auto& wsvs = workspace_variables();

  std::ostringstream os;

  os << "\n    Workspace& _ws [[maybe_unused]]";

  for (auto& t : wsm.out) {
    os << ",\n    std::optional<" << fix_type(wsvs.at(t).type) << ">& _" << t;
  }

  for (std::size_t i = 0; i < wsm.gout.size(); i++) {
    os << ",\n    std::optional<" << method_g_types(wsm.gout_type[i]) << ">& _"
       << wsm.gout[i];
  }

  for (auto& t : wsm.in) {
    if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
    os << ",\n    const std::optional<" << fix_type(wsvs.at(t).type) << ">& _" << t;
  }

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    os << ",\n    const std::optional<" << method_g_types(wsm.gin_type[i])
       << ">& _" << wsm.gin[i];
  }

  return os.str();
}

std::size_t count_any(const WorkspaceMethodInternalRecord& wsm) {
  return std::ranges::count(wsm.gout_type, "Any") +
         std::ranges::count(wsm.gin_type, "Any");
}

bool uses_variadic(const std::string& v) {
  return [](const std::string& s) { return s.find(',') != s.npos; }(v);
}

bool uses_variadic(const WorkspaceMethodInternalRecord& wsm) {
  return std::ranges::any_of(
             wsm.gout_type,
             [](const std::string& s) { return s.find(',') != s.npos; }) or
         std::ranges::any_of(wsm.gin_type, [](const std::string& s) {
           return s.find(',') != s.npos;
         });
}

std::string method_gout_selection(const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  for (std::size_t i = 0; i < wsm.gout.size(); i++) {
    if (wsm.gout_type[i] == "Any" or uses_variadic(wsm.gout_type[i])) {
      os << "      Wsv " << wsm.gout[i] << " = _" << wsm.gout[i] << " ? from(_"
         << wsm.gout[i]
         << R"(.value()) : throw std::runtime_error("Unknown output: \")"
         << wsm.gout[i] << "\\\"\");\n";
    } else {
      os << "      auto& " << wsm.gout[i] << " = select_gout<"
         << wsm.gout_type[i] << ">(_" << wsm.gout[i] << ", \"" << wsm.gout[i]
         << "\");\n";
    }
  }

  return os.str();
}

std::string method_gin_selection(const std::string& name,
                                 const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    const bool has_default = wsm.gin_value[i].has_value();

    if (wsm.gin_type[i] == "Any" or uses_variadic(wsm.gin_type[i])) {
      os << "      Wsv " << wsm.gin[i] << " = _" << wsm.gin[i]
         << " ? from(_" << wsm.gin[i]
         << R"(.value()) : throw std::runtime_error("Unknown input: \")"
         << wsm.gin[i] << "\\\"\");\n";
    } else {
      if (has_default) {
        os << "      static auto _" << wsm.gin[i]
           << "_default = []() { try { return workspace_methods().at(\"" << name
           << "\").defs.at(\"_" << wsm.gin[i] << "\").get<" << wsm.gin_type[i]
           << R"(>(); } catch(...) {throw std::runtime_error("DEV ERROR:\nFailed to initialize \")"
           << name << R"(\" default value \")" << wsm.gin[i]
           << "\\\"\"); } }();\n";
        os << "      auto& " << wsm.gin[i] << " = select_gin<"
           << wsm.gin_type[i] << ">(_" << wsm.gin[i] << ", _" << wsm.gin[i]
           << "_default);\n";
      } else {
        os << "      auto& " << wsm.gin[i] << " = select_gin<"
           << wsm.gin_type[i] << ">(_" << wsm.gin[i] << ", \"" << wsm.gin[i]
           << "\");\n";
      }
    }
  }

  return os.str();
}

std::string method_argument_selection(
    const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  static const auto& wsvs = workspace_variables();

  std::ostringstream os;

  for (auto& t : wsm.out) {
    os << "      auto& " << t << " = select_";
    if (std::ranges::any_of(wsm.in, Cmp::eq(t))) {
      os << "in";
    }
    os << "out<" << wsvs.at(t).type << ">(_" << t << ", _ws, \"" << t << "\");\n";
  }

  os << method_gout_selection(wsm);

  for (auto& t : wsm.in) {
    if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
    os << "      auto& " << t << " = select_in<"<<wsvs.at(t).type<<">(_"
       << t << ", _ws, \"" << t << "\");\n";
  }

  os << method_gin_selection(name, wsm);

  return os.str();
}

std::string unnamed_method_resolution_any(
    const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  std::size_t i_any = 0;
  for (std::size_t i = 0; i < wsm.gout.size(); i++) {
    if (wsm.gout_type[i] == "Any") {
      os << "      auto& _any" << ++i_any << " = " << wsm.gout[i] << ";\n";
    }
  }

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    if (wsm.gin_type[i] == "Any") {
      os << "      auto& _any" << ++i_any << " = " << wsm.gin[i] << ";\n";
    }
  }
  i_any = 1;

  os << "      std::visit([&](auto& _t) {"
        "\n        using _tT=std::remove_cvref_t<decltype(*_t)>;"
        "\n        "
     << name << "<_tT>(";

  bool any = false;
  if (wsm.pass_workspace) {
    os << "_ws";
    any = true;
  }

  for (auto& t : wsm.out) {
    if (any) os << ", ";
    any = true;
    os << t;
  }

  for (std::size_t i = 0; i < wsm.gout.size(); i++) {
    if (any) os << ", ";
    any = true;
    auto& t = wsm.gout[i];
    auto& tt = wsm.gout_type[i];
    if (tt == "Any") {
      if (i_any == 1)
        os << "*_t";
      else
        os << "_any" << i_any << ".get<_tT>()";
      i_any++;
    } else {
      os << t;
    }
  }

  for (auto& t : wsm.in) {
    if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
    if (any) os << ", ";
    any = true;
    os << t;
  }

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    if (any) os << ", ";
    any = true;
    auto& t = wsm.gin[i];
    auto& tt = wsm.gin_type[i];
    if (tt == "Any") {
      if (i_any == 1)
        os << "*_t";
      else
        os << "_any" << i_any << ".get<_tT>()";
      i_any++;
    } else {
      os << t;
    }
  }

  os << ");\n        }, _any1.value);\n";

  return os.str();
}

bool is_unique_variadic(const WorkspaceMethodInternalRecord& wsm) {
  for (auto& t : wsm.gout_type) {
    auto supergenerics = supergeneric_type(t);
    auto unique = unique_sorted(supergenerics);
    if (unique.size() != supergenerics.size()) return false;
    if (not std::is_permutation(
            unique.begin(), unique.end(), supergenerics.begin()))
      return false;
  }

  return true;
}

std::vector<std::string> unfix(std::vector<std::string> a) {
  for (auto& b: a) {
    if (b == "Numeric_") b = "std::shared_ptr<Numeric>";
    if (b == "Index_") b = "std::shared_ptr<Index>";
  }
  return a;
}

std::string unnamed_method_resolution_variadic(
    const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  const std::vector<std::vector<std::string>> supergenerics =
      supergeneric_types(wsm);
  if (is_unique_variadic(wsm)) {
    os << "      std::visit([&] (auto& _t) {"
          "\n        using _tT=std::remove_cvref_t<decltype(_t)>;";

    std::string var = "";
    os << "\n        if constexpr (false) {}";
    for (auto& generics : supergenerics) {
      std::size_t i_comma = 0;
      auto ungen = unfix(generics);

      os << "\n        else if constexpr (std::is_same_v<_tT, " << ungen[0]
         << ">) " << name << "(";

      bool first = true;
      for (auto& t : wsm.out) {
        if (not first) os << ", ";
        first = false;
        os << t;
      }

      for (std::size_t i = 0; i < wsm.gout.size(); i++) {
        if (not first) os << ", ";
        first = false;

        auto& t = wsm.gout_type[i];
        if (t.find(',') != t.npos) {
          if (i_comma == 0) {
            os << "*_t";
            var = wsm.gout[i];
          } else {
            os << "*std::get<" << ungen[i_comma] << ">(" << wsm.gout[i]
               << ".value)";
          }
          i_comma++;
        } else {
          os << wsm.gout[i];
        }
      }

      for (auto& t : wsm.in) {
        if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
        if (not first) os << ", ";
        first = false;
        os << t;
      }

      for (std::size_t i = 0; i < wsm.gin.size(); i++) {
        if (not first) os << ", ";
        first = false;

        auto& t = wsm.gin_type[i];
        if (t.find(',') != t.npos) {
          if (i_comma == 0) {
            os << "*_t";
            var = wsm.gin[i];
          } else {
            os << "*std::get<" << ungen[i_comma] << ">(" << wsm.gin[i]
               << ".value)";
          }
          i_comma++;
        } else {
          os << wsm.gin[i];
        }
      }

      os << ");";
    }
    os << "\n        else throw std::runtime_error(\"Bad combination in \\\"" << name << R"(\""); }, )" << var << ".value);\n";
  } else {
    os << "      ";
    for (auto& generics : supergenerics) {
      std::size_t i_comma = 0;
      auto ungen = unfix(generics);

      os << "if (";
      for (std::size_t i = 0; i < wsm.gout.size(); i++) {
        auto& t = wsm.gout_type[i];
        if (t.find(',') != t.npos) {
          if (i_comma > 0) os << " and ";
          os << "std::holds_alternative<";
          os << ungen[i_comma++];
          os << ">(" << wsm.gout[i] << ".value)";
        }
      }
      for (std::size_t i = 0; i < wsm.gin.size(); i++) {
        auto& t = wsm.gin_type[i];
        if (t.find(',') != t.npos) {
          if (i_comma > 0) os << " and ";
          os << "std::holds_alternative<";
          os << ungen[i_comma++];
          os << ">(" << wsm.gin[i] << ".value)";
        }
      }
      os << ") {\n        " << name << '(';

      bool first = true;
      for (auto& t : wsm.out) {
        if (not first) os << ", ";
        first = false;
        os << t;
      }
      i_comma = 0;
      for (std::size_t i = 0; i < wsm.gout.size(); i++) {
        if (not first) os << ", ";
        first = false;

        auto& t = wsm.gout_type[i];
        if (t.find(',') != t.npos) {
          os << "**std::get_if<" << ungen[i_comma] << ">(&" << wsm.gout[i]
             << ".value)";
          i_comma++;
        } else {
          os << wsm.gout[i];
        }
      }

      for (auto& t : wsm.in) {
        if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
        if (not first) os << ", ";
        first = false;
        os << t;
      }

      for (std::size_t i = 0; i < wsm.gin.size(); i++) {
        if (not first) os << ", ";
        first = false;

        auto& t = wsm.gin_type[i];
        if (t.find(',') != t.npos) {
          os << "**std::get_if<" << ungen[i_comma] << ">(&" << wsm.gin[i]
             << ".value)";
          i_comma++;
        } else {
          os << wsm.gin[i];
        }
      }

      os << ");\n      } else ";
    }
    os << "{\n";
    os << R"(        throw std::runtime_error("Cannot call method \")" << name
       << "\\\" with the given arguments,\");\n";
    os << "      }\n";
  }

  return os.str();
}

std::string unnamed_method_resolution_simple(
    const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  os << "      " << name << "(";

  bool any = false;
  if (wsm.pass_workspace) {
    os << "_ws";
    any = true;
  }

  for (auto& t : wsm.out) {
    if (any) os << ", ";
    any = true;
    os << t;
  }

  for (const auto& t : wsm.gout) {
    if (any) os << ", ";
    any = true;
    os << t;
  }

  for (auto& t : wsm.in) {
    if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
    if (any) os << ", ";
    any = true;
    os << t;
  }

  for (const auto& t : wsm.gin) {
    if (any) os << ", ";
    any = true;
    os << t;
  }

  return os.str() + ");\n";
}

std::string unnamed_method_resolution(
    const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  if (count_any(wsm) > 0) {
    os << '\n' << unnamed_method_resolution_any(name, wsm);
  } else if (uses_variadic(wsm)) {
    os << '\n' << unnamed_method_resolution_variadic(name, wsm);
  } else {
    os << '\n' << unnamed_method_resolution_simple(name, wsm);
  }

  return os.str();
}

std::string method_argument_documentation(
    const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  bool first = true;

  for (auto& t : wsm.out) {
    if (not first) os << ",\n    ";
    first = false;
    os << "py::arg_v(\"" << t << "\", std::nullopt, \"self." << t
       << "\").noconvert()";
  }

  for (auto& t : wsm.gout) {
    if (not first) os << ",\n    ";
    first = false;
    os << "py::arg_v(\"" << t << R"(", std::nullopt, "None").noconvert())";
  }

  for (auto& t : wsm.in) {
    if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
    if (not first) os << ",\n    ";
    first = false;
    os << "py::arg_v(\"" << t << "\", std::nullopt, \"self." << t << "\")";
  }

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    if (not first) os << ",\n    ";
    first = false;
    os << "py::arg_v(\"" << wsm.gin[i] << "\", std::nullopt, ";
    if (wsm.gin_value[i]) {
      os << "R\"-x-(" << to_defval_str(*wsm.gin_value[i], wsm.gin_type[i])
         << ")-x-\"";
    } else {
      os << R"("None")";
    }
    os << ")";
  }

  if (auto out = os.str(); out.size()) return out + ",\n    ";
  return "";
}

std::string method(const std::string& name,
                   const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  os << "  ws.def(\"" << name << "\",";
  os << "[](";
  os << method_arguments(wsm);
  os << ") -> void {\n";
  os << method_argument_selection(name, wsm);

  if (not wsm.pass_names)
    os << unnamed_method_resolution(name, wsm);
  else
    os << "\n      throw std::runtime_error(\"Cannot call method outside Agenda, see online documentation for workaround\");\n";

  os << "    },\n    " << method_argument_documentation(wsm) << "py::doc(R\""
     << method_docs(name) << "\"));\n\n";
  return os.str();
}

void methods(const std::string& fname) {
  const auto wsms = internal_workspace_methods();

  std::ofstream os(fname);

  os << R"--(#include <python_interface.h>

#include <m_copy.h>
#include <m_delete.h>
#include <m_general.h>
#include <m_ignore.h>
#include <m_xml.h>

namespace Python {
void py_auto_wsm(artsclass<Workspace>& ws [[maybe_unused]]) {
)--";

  for (auto& [name, wsv] : wsms) {
    os << method(name, wsv) << std::flush;
  }

  os << "}\n}  // namespace Python\n";
}

void groups(const std::string& fname) {
  const auto wsgs = internal_workspace_groups();

  std::ofstream hos(fname + ".h");

  hos << R"--(#include <auto_wsg.h>
#include <python_interface_value_type.h>

namespace Python {
using PyWsvValue = std::variant<
  )--";

  bool first = true;
  for (auto& [name, wsg] : wsgs) {
    if (not first) hos << ",\n    ";
    first = false;
    hos << "  ";
    if (name == "Index" or name == "Numeric") {
      hos << name << '_';
    } else {
      hos << "std::shared_ptr<" << name << '>';
    }
  }

  hos << R"--(>;

template <typename T>
concept PythonWorkspaceGroup =
       )--";

  first = true;
  for (auto& [name, wsg] : wsgs) {
    if (not first) hos << "\n    || ";
    first = false;
    hos << "std::is_same_v<T, " << name;
    if (name == "Index" or name == "Numeric") {
      hos << '_';
    }
    hos << '>';
  }

  hos << R"--(;

template <PythonWorkspaceGroup T>
Wsv from(const std::shared_ptr<T>& x) {
  Wsv out;
  if constexpr (WorkspaceGroup<T>)
    out.value = x;
  else
    out.value = x.val;
  return out;
}

template <WorkspaceGroup T>
Wsv from(const ValueHolder<T>& x) {
  Wsv out;
  out.value = x.val;
  return out;
}

template <typename T>
concept Wsvable = requires (T x) {
  { from(x) } -> std::same_as<Wsv>;
};

template <Wsvable ... Ts>
Wsv from(const std::variant<Ts...>& x) {
  return std::visit([](auto& arg) -> Wsv {
    return from(arg);
  }, x);
}

PyWsvValue from(const std::shared_ptr<Wsv>& x);
}  // namespace Python
)--";

  std::ofstream cos(fname + ".cpp");

  cos << R"--(#include <py_auto_wsg.h>

namespace Python {
PyWsvValue from(const std::shared_ptr<Wsv>& x) {
  return std::visit([](auto& arg) -> PyWsvValue {
    return arg;
  }, x -> value);
}
}  // namespace Python
)--";
}

int main() {
  groups("py_auto_wsg");
  variables("py_auto_wsv.cpp");
  methods("py_auto_wsm.cpp");
}
