#include <workspace.h>

#include <algorithm>
#include <iostream>
#include <ranges>
#include <sstream>
#include <string_view>

#include "compare.h"
#include "pydocs.h"

std::vector<std::string> errors;
#define ERRORAPPEND                \
  catch (std::exception & e) {     \
    errors.emplace_back(e.what()); \
  }

std::ofstream& select_ofstream(std::vector<std::ofstream>& ofs, int i) {
  return ofs[i % ofs.size()];
}

std::string using_pygroup() {
  std::ostringstream os;

  const auto& wsgs = internal_workspace_groups();

  for (auto& [name, wsg] : wsgs) {
    os << "  using py" << name << " [[maybe_unused]] = ";
    if (wsg.value_type) {
      os << "ValueHolder<" << name << ">";
    } else if (name == "Any") {
      os << "py::object";
    } else {
      os << name;
    }
    os << ";\n";
  }

  os << '\n';
  return os.str();
}

std::string method_arguments(const WorkspaceMethodInternalRecord& wsm) {
  const auto& wsvs = workspace_variables();

  std::ostringstream os;

  const auto generics = wsm.generic_overloads();

  os << "    Workspace& _ws [[maybe_unused]]";
  for (auto& v : wsm.out) {
    os << ",\n    py" << wsvs.at(v).type << "* const _" << v;
  }

  // NB: This is untested code because we do not have any methods with generic outputs at time of writing
  for (Size i = 0; i < wsm.gout.size(); i++) {
    const auto& v = wsm.gout.at(i);

    if (generics[i].size() == 1) {
      const auto& g = wsm.gout_type.at(i);
      os << ",\n    py" << g << "* const _" << v;
    } else {
      const auto& v          = wsm.gin.at(i);
      std::string_view comma = "";
      os << ",\n    std::variant<";
      for (auto& arg : generics[i]) {
        os << std::exchange(comma, ", ") << "std::shared_ptr<py" << arg << ">";
      }
      os << "> * const _" << v;
    }
  }

  const auto not_out = [&wsm](const auto& v) {
    return std::ranges::none_of(wsm.out, Cmp::eq(v));
  };

  for (auto& v : wsm.in | std::ranges::views::filter(not_out)) {
    os << ",\n    const py" << wsvs.at(v).type << "* const _" << v;
  }

  const auto is_gout = [&wsm](const auto& v) {
    return std::ranges::any_of(wsm.gout, Cmp::eq(v));
  };

  // NB: This is partly untested code because we do not have any methods with generic input+output at time of writing
  Size index = wsm.gout.size();
  for (Size i = 0; i < wsm.gin.size(); i++) {
    const auto& v = wsm.gin.at(i);

    if (is_gout(v)) continue;

    if (generics[index].size() == 1) {
      const auto& g = wsm.gin_type.at(i);
      os << ",\n    const py" << g << "* const _" << v;
    } else {
      const auto& v          = wsm.gin.at(i);
      std::string_view comma = "";
      os << ",\n    const std::variant<";
      for (auto& arg : generics[index]) {
        os << std::exchange(comma, ", ") << "std::shared_ptr<py" << arg << ">";
      }
      os << "> * const _" << v;
    }

    index++;
  }

  return os.str();
}

bool uses_variadic(const std::string& v) {
  return [](const std::string& s) { return s.find(',') != s.npos; }(v);
}

bool uses_variadic(const WorkspaceMethodInternalRecord& wsm) {
  return std::ranges::any_of(
             wsm.gout_type,
             [](const std::string& s) { return s.find(',') != s.npos; }) or
         std::ranges::any_of(wsm.gin_type, [](const std::string& s) {
           return uses_variadic(s);
         });
}

std::string method_gout_selection(const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  for (std::size_t i = 0; i < wsm.gout.size(); i++) {
    if (wsm.gout_type[i] == "Any" or uses_variadic(wsm.gout_type[i])) {
      os << "        Wsv " << wsm.gout[i] << " = _" << wsm.gout[i]
         << " ? from(_" << wsm.gout[i]
         << R"() : throw std::runtime_error("Unknown output: \")" << wsm.gout[i]
         << "\\\"\");\n";
    } else {
      os << "        auto& " << wsm.gout[i] << " = select_gout<"
         << wsm.gout_type[i] << ">(_" << wsm.gout[i] << ", _ws, \""
         << wsm.gout[i] << "\");\n";
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
      os << "        Wsv " << wsm.gin[i] << " = _" << wsm.gin[i] << " ? from(_"
         << wsm.gin[i] << R"() : throw std::runtime_error("Unknown input: \")"
         << wsm.gin[i] << "\\\"\");\n";
    } else {
      if (has_default) {
        os << "        static const auto _" << wsm.gin[i]
           << "_default = []() { try { return workspace_methods().at(\"" << name
           << "\").defs.at(\"_" << wsm.gin[i] << "\").get<" << wsm.gin_type[i]
           << R"(>(); } catch(...) {throw std::runtime_error("DEV ERROR:\nFailed to initialize \")"
           << name << R"(\" default value \")" << wsm.gin[i]
           << "\\\"\"); } }();\n";
        os << "        auto& " << wsm.gin[i] << " = select_gin<"
           << wsm.gin_type[i] << ">(_" << wsm.gin[i] << ", _" << wsm.gin[i]
           << "_default);\n";
      } else {
        os << "        auto& " << wsm.gin[i] << " = select_gin<"
           << wsm.gin_type[i] << ">(_" << wsm.gin[i] << ", \"" << wsm.gin[i]
           << "\");\n";
      }
    }
  }

  return os.str();
}

std::string method_argument_selection(
    const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  const auto& wsvs = workspace_variables();

  std::ostringstream os;

  for (auto& t : wsm.out) {
    os << "        auto& " << t << " = select_";
    if (std::ranges::any_of(wsm.in, Cmp::eq(t))) {
      os << "in";
    }
    os << "out<" << wsvs.at(t).type << ">(_" << t << ", _ws, \"" << t
       << "\");\n";
  }

  os << method_gout_selection(wsm);

  for (auto& t : wsm.in) {
    if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
    os << "        auto& " << t << " = select_in<" << wsvs.at(t).type << ">(_"
       << t << ", _ws, \"" << t << "\");\n";
  }

  os << method_gin_selection(name, wsm);

  return os.str();
}

std::string method_resolution_any(const std::string& name,
                                  const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  std::size_t i_any = 0;
  for (std::size_t i = 0; i < wsm.gout.size(); i++) {
    if (wsm.gout_type[i] == "Any") {
      os << "        auto& _any" << ++i_any << " = " << wsm.gout[i] << ";\n";
    }
  }

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    if (wsm.gin_type[i] == "Any") {
      os << "        auto& _any" << ++i_any << " = " << wsm.gin[i] << ";\n";
    }
  }
  i_any = 1;

  os << "        try {\n";
  os << "          std::visit([&](auto& _t) {"
        "\n           using _tT=std::remove_cvref_t<decltype(*_t)>;"
        "\n            "
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
    any      = true;
    auto& t  = wsm.gout[i];
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
    any      = true;
    auto& t  = wsm.gin[i];
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

  os << ");\n          }, _any1.value());\n";
  os << "        } catch (std::bad_variant_access&) {\n";

  const auto* const spaces = "                                              ";
  os << R"(          throw std::runtime_error(var_string("Derived variadic arguments:\n",)";

  for (std::size_t i = 0; i < wsm.gout.size(); i++) {
    auto& t = wsm.gout_type[i];
    auto& v = wsm.gout[i];
    if (uses_variadic(t) or t == "Any") {
      os << '\n'
         << spaces << "\"" << v << " : \", " << v << ".type_name(), '\\n',";
    }
  }

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    auto& t = wsm.gin_type[i];
    auto& v = wsm.gin[i];
    if (uses_variadic(t) or t == "Any") {
      os << '\n'
         << spaces << "\"" << v << " : \", " << v << ".type_name(), '\\n',";
    }
  }

  os << '\n'
     << spaces
     << "\"If the signature above give different types for any of these variable(s),\\n\",\n"
     << spaces
     << "\"then you need to specify the type of the variable(s) explicitly.\\n\",\n"
     << spaces
     << "\"Otherwise, the desired type combination is not available for this method.\"";
  os << "));\n";

  os << "        } catch (std::exception&) { throw; }\n";

  return os.str();
}

std::vector<std::string> supergeneric_type(const std::string& t) {
  std::vector<std::string> types;
  std::string type = t;
  auto ptr         = type.find(',');
  while (ptr != type.npos) {
    types.push_back(type.substr(0, ptr));
    type = type.substr(ptr + 1);
    ptr  = type.find(',');
  }
  types.push_back(type);

  for (auto& T : types) {
    while (T.back() == ' ') T.pop_back();
    while (T.front() == ' ') T.erase(0, 1);
  }

  return types;
}

std::vector<std::string> unique_sorted(const std::vector<std::string>& v) {
  std::vector<std::string> out = v;
  std::ranges::sort(out);
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

bool is_unique_variadic(const WorkspaceMethodInternalRecord& wsm) {
  for (auto& t : wsm.gout_type) {
    auto supergenerics = supergeneric_type(t);
    auto unique        = unique_sorted(supergenerics);
    if (unique.size() != supergenerics.size()) return false;
    if (not std::is_permutation(
            unique.begin(), unique.end(), supergenerics.begin()))
      return false;
  }

  return true;
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

std::string method_resolution_variadic(
    const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  const std::vector<std::vector<std::string>> supergenerics =
      supergeneric_types(wsm);
  if (is_unique_variadic(wsm)) {
    os << "        std::visit([&] (auto& _t) {"
          "\n          using _tT=std::remove_cvref_t<decltype(*_t)>;";

    std::string var = "";
    os << "\n          if constexpr (false) {}";
    for (auto& generics : supergenerics) {
      std::size_t i_comma = 0;
      auto ungen          = generics;

      os << "\n          else if constexpr (std::is_same_v<_tT, " << ungen[0]
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
            os << wsm.gout[i] << ".get<" << ungen[i_comma] << ">()";
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
            os << wsm.gin[i] << ".get<" << ungen[i_comma] << ">()";
          }
          i_comma++;
        } else {
          os << wsm.gin[i];
        }
      }

      os << ");";
    }

    const auto* const spaces =
        "                                                    ";
    os << "\n          else throw std::runtime_error(var_string(\"Derived variadic arguments:\\n\",";

    for (std::size_t i = 0; i < wsm.gout.size(); i++) {
      auto& t = wsm.gout_type[i];
      auto& v = wsm.gout[i];
      if (uses_variadic(t) or t == "Any") {
        os << '\n'
           << spaces << "\"" << v << " : \", " << v << ".type_name(), '\\n',";
      }
    }

    for (std::size_t i = 0; i < wsm.gin.size(); i++) {
      auto& t = wsm.gin_type[i];
      auto& v = wsm.gin[i];
      if (uses_variadic(t) or t == "Any") {
        os << '\n'
           << spaces << "\"" << v << " : \", " << v << ".type_name(), '\\n',";
      }
    }

    os << '\n'
       << spaces
       << "\"If the signature above give different types for any of these variable(s),\\n\",\n"
       << spaces
       << "\"then you need to specify the type of the variable(s) explicitly.\\n\",\n"
       << spaces
       << "\"Otherwise, the desired type combination is not available for this method.\"";
    os << "));\n"
          "        }, "
       << var << ".value());\n";
  } else {
    os << "      ";
    for (auto& generics : supergenerics) {
      std::size_t i_comma = 0;
      auto ungen          = generics;

      os << "if (";
      for (std::size_t i = 0; i < wsm.gout.size(); i++) {
        auto& t = wsm.gout_type[i];
        if (t.find(',') != t.npos) {
          if (i_comma > 0) os << " and ";
          os << wsm.gout[i] << ".holds<" << ungen[i_comma++] << ">()";
        }
      }
      for (std::size_t i = 0; i < wsm.gin.size(); i++) {
        auto& t = wsm.gin_type[i];
        if (t.find(',') != t.npos) {
          if (i_comma > 0) os << " and ";
          os << wsm.gin[i] << ".holds<" << ungen[i_comma++] << ">()";
        }
      }
      os << ") {\n          " << name << '(';

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
          os << wsm.gout[i] << ".get_unsafe<" << ungen[i_comma] << ">()";
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
          os << wsm.gin[i] << ".get_unsafe<" << ungen[i_comma] << ">()";
          i_comma++;
        } else {
          os << wsm.gin[i];
        }
      }

      os << ");\n        } else ";
    }
    os << "{\n";

    const auto* const spaces =
        "                                               ";
    os << "\n          throw std::runtime_error(var_string(\"Derived variadic arguments:\\n\",";

    for (std::size_t i = 0; i < wsm.gout.size(); i++) {
      auto& t = wsm.gout_type[i];
      auto& v = wsm.gout[i];
      if (uses_variadic(t) or t == "Any") {
        os << '\n'
           << spaces << "\"" << v << " : \", " << v << ".type_name(), '\\n',";
      }
    }

    for (std::size_t i = 0; i < wsm.gin.size(); i++) {
      auto& t = wsm.gin_type[i];
      auto& v = wsm.gin[i];
      if (uses_variadic(t) or t == "Any") {
        os << '\n'
           << spaces << "\"" << v << " : \", " << v << ".type_name(), '\\n',";
      }
    }

    os << '\n'
       << spaces
       << "\"If the signature above give different types for any of these variable(s),\\n\",\n"
       << spaces
       << "\"then you need to specify the type of the variable(s) explicitly.\\n\",\n"
       << spaces
       << "\"Otherwise, the desired type combination is not available for this method.\"";
    os << "));\n";
    os << "        }\n";
  }

  return os.str();
}

std::string method_resolution_simple(const std::string& name,
                                     const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  os << "        " << name << "(";

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

std::size_t count_any(const WorkspaceMethodInternalRecord& wsm) {
  return std::ranges::count(wsm.gout_type, "Any") +
         std::ranges::count(wsm.gin_type, "Any");
}

std::string method_resolution(const std::string& name,
                              const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  if (count_any(wsm) > 0) {
    os << '\n' << method_resolution_any(name, wsm);
  } else if (uses_variadic(wsm)) {
    os << '\n' << method_resolution_variadic(name, wsm);
  } else {
    os << '\n' << method_resolution_simple(name, wsm);
  }

  return os.str();
}

std::string method_error(const std::string& name,
                         const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;
  os << "        throw std::runtime_error(\n          var_string(";
  os << R"("Cannot execute method:\n\n",
              ")";
  os << name << '(';

  const auto spaces = std::string(name.size() + 1, ' ');

  bool first = true;
  for (auto& t : wsm.out) {
    if (not first) os << ",\\n" << spaces;
    first = false;
    os << t << " : \",\n              "
       << "_" << t << " ? type(_" << t << R"() : std::string(R"-x-()"
       << workspace_variables().at(t).type;
    os << R"(, defaults to self.)" << t << R"()-x-"),
              ")";
  }

  for (auto& t : wsm.gout) {
    if (not first) os << ",\\n" << spaces;
    first = false;
    os << t << " : \",\n              "
       << "_" << t << " ? type(_" << t << R"() : std::string("NoneType"),
              ")";
  }

  for (auto& t : wsm.in) {
    if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
    if (not first) os << ",\\n" << spaces;
    first = false;
    os << t << " : \",\n              "
       << "_" << t << " ? type(_" << t << R"() : std::string(R"-x-()"
       << workspace_variables().at(t).type;
    os << R"(, defaults to self.)" << t << R"()-x-"),
              ")";
  }

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    auto& t = wsm.gin[i];
    auto& v = wsm.gin_value[i];
    auto& g = wsm.gin_type[i];
    if (not first) os << ",\\n" << spaces;
    first = false;
    os << t << " : \",\n              "
       << "_" << t << " ? type(_" << t << R"() : std::string(R"-x-()";
    if (v) {
      os << g << ", defaults to " << to_defval_str(*v);
    } else {
      os << "NoneType";
    }
    os << R"()-x-"),
              ")";
  }

  os << ")\", ";
  os << "\"\\n\\n\", e.what()));\n";
  return os.str();
}

std::string method_argument_documentation(
    const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  bool first = true;

  for (const auto& t : wsm.out) {
    if (not first) os << ",\n    ";
    first = false;
    os << "\"" << t << "\"_a.noconvert().none() = py::none()";
  }

  for (const auto& t : wsm.gout) {
    if (not first) os << ",\n    ";
    first = false;
    os << "\"" << t << "\"_a.noconvert().none() = py::none()";
  }

  for (const auto& t : wsm.in) {
    if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
    if (not first) os << ",\n    ";
    first = false;
    os << "\"" << t << "\"_a.none() = py::none()";
  }

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    if (not first) os << ",\n    ";
    first = false;
    if (wsm.gin_value[i]) {
      os << "\"" << wsm.gin[i] << "\"_a.none() = py::none()";
    } else {
      os << "\"" << wsm.gin[i] << "\"_a.none() = py::none()";
    }
  }

  if (auto out = os.str(); out.size()) return out + ",\n    ";
  return "";
}

std::string method(const std::string& name,
                   const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  os << "  ws.def(\"" << name << "\",";
  os << "[](\n";
  os << method_arguments(wsm);
  os << ") -> void {\n      try {\n";
  os << method_argument_selection(name, wsm);
  os << method_resolution(name, wsm);

  os << "      } catch (std::exception& e) {\n";
  os << method_error(name, wsm);
  os << "      }\n";
  os << "    },\n    " << method_argument_documentation(wsm) << "R\""
     << method_docs(name) << "\",\n";
  os << "    py::call_guard<py::gil_scoped_release>());\n\n";
  return os.str();
}

void methods(int nfiles) {
  const auto& wsms = internal_workspace_methods();
  const auto wsas  = internal_workspace_agendas();

  std::vector<std::ofstream> ofs(nfiles);
  for (int i = 0; i < nfiles; i++) {
    ofs[i].open(std::format("py_auto_wsm_{}.cpp", i));
  }

  for (int i = 0; i < nfiles; i++) {
    select_ofstream(ofs, i) << R"--(#include <python_interface.h>

#include <m_ignore.h>
#include <m_xml.h>
#include <workspace.h>

#include <nanobind/stl/variant.h>
#include <nanobind/stl/shared_ptr.h>

namespace Python {
void py_auto_wsm_)--" << i << "(py::class_<Workspace>& ws [[maybe_unused]]) {\n"
                            << using_pygroup();
  }

  int ifile = 0;
  for (auto& [name, wsv] : wsms) {
    try {
      select_ofstream(ofs, ifile++) << method(name, wsv) << std::flush;
    }
    ERRORAPPEND;
  }

  for (int i = 0; i < nfiles; i++) {
    select_ofstream(ofs, i) << "}\n}  // namespace Python\n";
  }
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " output_file_count[int]" << std::endl;
    return EXIT_FAILURE;
  }

  const int num_methods = std::stoi(argv[1]);

  methods(num_methods);

  if (errors.size()) {
    std::cerr << "Errors (" << errors.size() << "):\n";
    for (auto& e : errors) {
      std::cerr << e << '\n';
    }
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
