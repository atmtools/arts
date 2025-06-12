#include <workspace.h>

#include <algorithm>
#include <iostream>
#include <ranges>
#include <sstream>
#include <string_view>

#include "compare.h"
#include "pydocs.h"
#include "workspace_groups.h"

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
    if (wsm.gout_type[i] == "Any") {
      std::println(
          os,
          R"-x-(        Wsv {0} = _{0} ? from(_{0}) : throw std::runtime_error(R"(Unknown output: "{0}")");)-x-",
          wsm.gout[i]);
    } else if (uses_variadic(wsm.gout_type[i])) {
      std::println(
          os,
          R"-x-(        auto& {0} = _{0} ? *_{0} : throw std::runtime_error(R"(Unknown output: "{0}")");)-x-",
          wsm.gout[i],
          wsm.gout_type[i]);
    } else {
      std::println(os,
                   R"(        {1}& {0} = select_gout<{1}>(_{0}, _ws, "{0}");)",
                   wsm.gout[i],
                   wsm.gout_type[i]);
    }
  }

  return os.str();
}

std::string method_gin_selection(const std::string& name,
                                 const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    const bool has_default = wsm.gin_value[i].has_value();

    if (wsm.gin_type[i] == "Any") {
      std::println(
          os,
          R"-x-(        const Wsv {0} = _{0} ? from(_{0}) : throw std::runtime_error(R"(Unknown input: "{0}")");)-x-",
          wsm.gin[i]);
    } else if (uses_variadic(wsm.gin_type[i])) {
      std::println(
          os,
          R"-x-(        const auto& {0} = _{0} ? *_{0} : throw std::runtime_error(R"(Unknown input: "{0}")");)-x-",
          wsm.gin[i]);
    } else {
      if (has_default) {
        std::println(os,
                     R"x(        static const {2} _{0}_default = []() -> {2} {{
            try {{
              return workspace_methods().at("{1}").defs.at("_{0}").get<{2}>();
            }} catch(...) {{
              throw std::runtime_error(R"(DEV ERROR:
Failed to initialize "{1}" default value "_{0}")");
            }}
          }}();
        const {2}& {0} = select_gin<{2}>(_{0}, _{0}_default);)x",
                     wsm.gin[i],
                     name,
                     wsm.gin_type[i]);
      } else {
        std::println(
            os,
            R"x(        const {1}& {0} = select_gin<{1}>(_{0}, "{0}");)x",
            wsm.gin[i],
            wsm.gin_type[i]);
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
    std::println(
        os,
        R"x(        {2}& {0} = select_{1}out<{2}>(_{0}, _ws, "{0}");)x",
        t,
        stdr::any_of(wsm.in, Cmp::eq(t)) ? "in"sv : ""sv,
        wsvs.at(t).type);
  }

  os << method_gout_selection(wsm);

  for (auto& t : wsm.in) {
    if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
    std::println(
        os,
        R"x(        const {1}& {0} = select_in<{1}>(_{0}, _ws, "{0}");)x",
        t,
        wsvs.at(t).type);
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
      std::println(os, "        auto& _any{} = {};", ++i_any, wsm.gout[i]);
    }
  }

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    if (wsm.gin_type[i] == "Any") {
      std::println(os, "        auto& _any{} = {};", ++i_any, wsm.gin[i]);
    }
  }

  std::println(os, "        switch(_any1.value_index()) {{");
  for (auto& group : internal_workspace_groups() | stdv::keys) {
    i_any = 1;

    std::print(
        os,
        R"(          case WorkspaceGroupInfo<{0}>::index: return {1}<{0}>({2})",
        group,
        name,
        wsm.pass_workspace ? "_ws"sv : ""sv);

    bool any = wsm.pass_workspace;

    for (auto& t : wsm.out) {
      std::print(os, "{0}{1}", any ? ", "sv : ""sv, t);
      any = true;
    }

    for (std::size_t i = 0; i < wsm.gout.size(); i++) {
      auto& t = wsm.gout[i];
      if (wsm.gout_type[i] == "Any") {
        std::print(os,
                   "{3}_any{0}.get{1}<{2}>()",
                   i_any,
                   i_any == 1 ? "_unsafe"sv : ""sv,
                   group,
                   any ? ", "sv : ""sv);
        i_any++;
      } else {
        std::print(os, "{0}{1}", any ? ", "sv : ""sv, t);
      }
      any = true;
    }

    for (auto& t : wsm.in) {
      if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
      std::print(os, "{0}{1}", any ? ", "sv : ""sv, t);
      any = true;
    }

    for (std::size_t i = 0; i < wsm.gin.size(); i++) {
      auto& t = wsm.gin[i];
      if (wsm.gin_type[i] == "Any") {
        std::print(os,
                   "{3}_any{0}.get{1}<{2}>()",
                   i_any,
                   i_any == 1 ? "_unsafe"sv : ""sv,
                   group,
                   any ? ", "sv : ""sv);
        i_any++;
      } else {
        std::print(os, "{0}{1}", any ? ", "sv : ""sv, t);
      }
      any = true;
    }

    std::println(os, ");");
  }
  std::println(
      os,
      "        }}\n        throw std::runtime_error(\"Cannot understand input\");");

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
    std::string var = "";
    for (auto& generics : supergenerics) {
      const auto& ungen = generics;

      std::string test{"true"};
      std::vector<std::string> vars{};

      for (std::size_t i = 0; i < wsm.gout.size(); i++) {
        auto& t = wsm.gout_type[i];
        if (uses_variadic(t)) {
          const auto& mytype = ungen[vars.size()];
          std::format_to(
              std::back_inserter(test),
              " and std::holds_alternative<{0}>(std::shared_ptr<py{1}>)",
              mytype,
              wsm.gout[i]);
          vars.emplace_back(
              std::format("**std::get_if<py{0}>({1})", mytype, wsm.gout[i]));
        }
      }

      for (std::size_t i = 0; i < wsm.gin.size(); i++) {
        auto& t = wsm.gin_type[i];
        if (uses_variadic(t)) {
          const auto& mytype = ungen[vars.size()];
          std::format_to(
              std::back_inserter(test),
              " and std::holds_alternative<std::shared_ptr<py{0}>>({1})",
              mytype,
              wsm.gin[i]);
          vars.emplace_back(
              std::format("**std::get_if<std::shared_ptr<py{0}>>(&{1})",
                          mytype,
                          wsm.gin[i]));
        }
      }

      stdr::reverse(vars);
      std::print(os,
                 "        if ({0})\n          return {1}({2}",
                 test,
                 name,
                 wsm.pass_workspace ? "_ws"sv : ""sv);

      bool any = wsm.pass_workspace;

      for (auto& t : wsm.out) {
        std::print(os, "{0}{1}", any ? ", "sv : ""sv, t);
        any = true;
      }

      for (std::size_t i = 0; i < wsm.gout.size(); i++) {
        if (uses_variadic(wsm.gout_type[i])) {
          std::print(os, "{0}{1}", any ? ", "sv : ""sv, vars.back());
          vars.pop_back();
        } else {
          std::print(os, "{0}{1}", any ? ", "sv : ""sv, wsm.gout[i]);
        }
        any = true;
      }

      for (auto& t : wsm.in) {
        if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
        std::print(os, "{0}{1}", any ? ", "sv : ""sv, t);
        any = true;
      }

      for (std::size_t i = 0; i < wsm.gin.size(); i++) {
        if (uses_variadic(wsm.gin_type[i])) {
          std::print(os, "{0}{1}", any ? ", "sv : ""sv, vars.back());
          vars.pop_back();
        } else {
          std::print(os, "{0}{1}", any ? ", "sv : ""sv, wsm.gin[i]);
        }
        any = true;
      }

      os << ");\n";
    }

    os << R"(        throw std::runtime_error("Type mismatch.");)";
  } else {
    throw std::runtime_error("Not implemented");
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

Size max_varlen(const WorkspaceMethodInternalRecord& wsm) {
  using std::max;
  Size n = 0;
  for (auto& t : wsm.out) n = max(n, t.size());
  for (auto& t : wsm.in) n = max(n, t.size());
  for (auto& t : wsm.gin) n = max(n, t.size());
  for (auto& t : wsm.gout) n = max(n, t.size());
  return n;
}

std::string method_error(const std::string& name,
                         const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;
  os << std::format(R"(        throw std::runtime_error(
          std::format(R"-WSM-(Cannot execute method:

  {}()",
                    name);

  const auto spaces      = std::string(name.size() + 1 + 2, ' ');
  const Size largest_var = max_varlen(wsm);
  std::vector<std::string> arg;

  bool first = true;
  for (auto& t : wsm.out) {
    if (not first) os << ",\n" << spaces;
    first = false;
    std::print(
        os, R"({1} {0}: {{}})", std::string(largest_var - t.size(), ' '), t);
    arg.push_back(std::format(R"(_{0} ? "User-provided {1}"sv : "self.{0}"sv)",
                              t,
                              workspace_variables().at(t).type));
  }

  for (std::size_t i = 0; i < wsm.gout.size(); i++) {
    auto& t  = wsm.gout[i];
    auto& tt = wsm.gout_type[i];
    if (not first) os << ",\n" << spaces;
    first = false;
    std::print(
        os, R"({1} {0}: {{}})", std::string(largest_var - t.size(), ' '), t);
    if (uses_variadic(tt) or tt == "Any") {
      arg.push_back(std::format(
          R"(_{0} ? std::format("User-provided {{}}", type(_{0})) : std::string("None"))",
          t));
    } else {
      arg.push_back(
          std::format(R"(_{0} ? "User-provided {1}"sv : "self.{0}"sv)", t, tt));
    }
  }

  for (auto& t : wsm.in) {
    if (std::ranges::any_of(wsm.out, Cmp::eq(t))) continue;
    if (not first) os << ",\n" << spaces;
    first = false;
    std::print(
        os, R"({1} {0}: {{}})", std::string(largest_var - t.size(), ' '), t);
    arg.push_back(std::format(R"(_{0} ? "User-provided {1}"sv : "self.{0}"sv)",
                              t,
                              workspace_variables().at(t).type));
  }

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    auto& t = wsm.gin[i];
    auto& v = wsm.gin_value[i];
    if (not first) os << ",\n" << spaces;
    first = false;
    std::print(
        os, R"({1} {0}: {{}})", std::string(largest_var - t.size(), ' '), t);
    arg.push_back(std::format(
        R"(_{0} ? std::format("User-provided {{}}", type(_{0})) : std::string(R"-WSMVAR-({1})-WSMVAR-"))",
        t,
        v ? std::format("{}", to_defval_str(*v)) : "None"));
  }

  os << ")\n\nMethod reports the following error(s):\n{})-WSM-\", ";
  for (auto& a : arg) {
    os << std::format("{},\n           ", a);
  }
  os << "std::string_view(e.what())));\n";
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
  return std::format(
      R"-x-(  ws.def("{0}",[]({1}) -> void {{
    try {{
{2}{3}
    }} catch (std::exception& e) {{
{4}      }}
    }},
    {5}
{6},
    py::call_guard<py::gil_scoped_release>());

)-x-",
      name,
      method_arguments(wsm),
      method_argument_selection(name, wsm),
      method_resolution(name, wsm),
      method_error(name, wsm),
      method_argument_documentation(wsm),
      method_docs(name));
}

void methods(int nfiles) {
  const auto& wsms = internal_workspace_methods();

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
    std::cerr << "Usage: " << argv[0] << " output_file_count[int]" << '\n';
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
