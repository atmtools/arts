#include <workspace.h>
#include <workspace_dimensions.h>

#include <algorithm>
#include <iostream>
#include <ranges>
#include <sstream>
#include <string_view>

#include "compare.h"
#include "pydocs.h"
#include "workspace_groups.h"

namespace {
std::vector<std::string> errors;

#define ERRORAPPEND                \
  catch (std::exception & e) {     \
    errors.emplace_back(e.what()); \
  }

std::ofstream& select_ofstream(std::vector<std::ofstream>& ofs, int i) { return ofs[i % ofs.size()]; }

std::string using_pygroup() {
  std::ostringstream os;

  const auto& wsgs = internal_workspace_groups();

  for (auto& [name, wsg] : wsgs) {
    os << "  using py" << name << " [[maybe_unused]] = ";
    if (wsg.value_type) {
      os << "ValueHolder<" << name << ">";
    } else if (name == "Any") {
      os << "Wsv";
    } else {
      os << name;
    }
    os << ";\n";
  }

  os << '\n';
  return os.str();
}

std::string method_arguments(const WorkspaceMethodInternalRecord& wsm) {
  const auto&        wsvs = workspace_variables();
  std::ostringstream os;
  os << "    Workspace& _ws [[maybe_unused]]";
  for (const auto& name : wsm.out) os << ",\n    py" << wsvs.at(name).type << "* const _" << name;
  for (Size i = 0; i < wsm.gout.size(); ++i) {
    const auto& type = wsm.gout_type[i];
    os << ",\n    " << (type == "Any" or type.contains(',') ? "const py::object" : "py" + type) << "* const _"
       << wsm.gout[i];
  }
  for (const auto& name : wsm.in) {
    if (stdr::find(wsm.out, name) != wsm.out.end()) continue;
    os << ",\n    const py" << wsvs.at(name).type << "* const _" << name;
  }
  for (Size i = 0; i < wsm.gin.size(); ++i) {
    if (stdr::find(wsm.gout, wsm.gin[i]) != wsm.gout.end()) continue;
    const auto& type = wsm.gin_type[i];
    os << ",\n    const " << (type == "Any" or type.contains(',') ? "py::object" : "py" + type) << "* const _"
       << wsm.gin[i];
  }
  return os.str();
}

std::string generic_selection(const std::string& name,
                              const std::string& declaration,
                              bool               output,
                              const std::string& default_value = "") {
  const auto  type = WorkspaceMethodInternalRecord::generic_type(declaration, output);
  std::string allowed;
  if (declaration != "Any") {
    for (auto item : split(declaration, ",")) {
      trim(item);
      if (not allowed.empty()) allowed += ", ";
      allowed += "\"" + item + "\"";
    }
  }
  auto conversion = std::format("from_allowed(_{0}, {{{1}}}, {2})", name, allowed, output ? "true" : "false");
  if (not default_value.empty())
    conversion = std::format("(_{0} and not _{0}->is_none()) ? {1} : {2}", name, conversion, default_value);
  return std::format(
      "        auto {0}_owner = {1};\n        auto {0} = {2}::from({0}_owner);\n", name, conversion, type);
}

bool uses_variadic(const std::string& v) {
  return [](const std::string& s) { return s.find(',') != s.npos; }(v);
}

std::string method_gout_selection(const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  for (std::size_t i = 0; i < wsm.gout.size(); i++) {
    if (wsm.gout_type[i] == "Any" or uses_variadic(wsm.gout_type[i])) {
      os << generic_selection(wsm.gout[i], wsm.gout_type[i], true);
    } else {
      std::println(os, R"(        {1}& {0} = select_gout<{1}>(_{0}, _ws, "{0}");)", wsm.gout[i], wsm.gout_type[i]);
    }
  }

  return os.str();
}

std::string method_gin_selection(const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    const bool has_default = wsm.gin_value[i].has_value();

    if (stdr::find(wsm.gout, wsm.gin[i]) != wsm.gout.end()) continue;
    if (wsm.gin_type[i] == "Any" or uses_variadic(wsm.gin_type[i])) {
      const auto fallback =
          has_default ? std::format("workspace_methods().at(\"{}\").defs.at(\"_{}\")", name, wsm.gin[i]) : "";
      os << generic_selection(wsm.gin[i], wsm.gin_type[i], false, fallback);
    } else {
      if (has_default) {
        std::println(os,
                     R"x(        static const {2} _{0}_default = [] -> {2} {{
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
        std::println(os, R"x(        const {1}& {0} = select_gin<{1}>(_{0}, "{0}");)x", wsm.gin[i], wsm.gin_type[i]);
      }
    }
  }

  return os.str();
}

std::string method_argument_selection(const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  const auto& wsvs = workspace_variables();

  std::ostringstream os;

  for (auto& t : wsm.out) {
    std::println(os,
                 R"x(        {2}& {0} = select_{1}out<{2}>(_{0}, _ws, "{0}");)x",
                 t,
                 stdr::any_of(wsm.in, Cmp::eq(t)) ? "in"sv : ""sv,
                 wsvs.at(t).type);
  }

  os << method_gout_selection(wsm);

  for (auto& t : wsm.in) {
    if (stdr::any_of(wsm.out, Cmp::eq(t))) continue;
    std::println(os, R"x(        const {1}& {0} = select_in<{1}>(_{0}, _ws, "{0}");)x", t, wsvs.at(t).type);
  }

  os << method_gin_selection(name, wsm);

  return os.str();
}

std::string method_output_checks(const WorkspaceMethodInternalRecord& wsm) {
  std::string out;
  out += size_check_code(method_output_size_checks(wsm), "        ");
  return out;
}

std::string method_resolution_simple(const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  // Anything the method created has to be verified before returning it, so the
  // call cannot be the return statement when there is something to verify
  const auto checks    = wsm.return_type == "void" ? method_output_checks(wsm) : std::string{};
  const auto returning = checks.empty();

  os << (returning ? "        return " : "        ") << name << "(";

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
    if (stdr::any_of(wsm.out, Cmp::eq(t))) continue;
    if (any) os << ", ";
    any = true;
    os << t;
  }

  for (const auto& t : wsm.gin) {
    if (stdr::find(wsm.gout, t) != wsm.gout.end()) continue;
    if (any) os << ", ";
    any = true;
    os << t;
  }

  os << ");\n";

  if (not returning) os << checks << "        return;\n";

  return os.str();
}

std::string method_resolution(const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  return method_resolution_simple(name, wsm);
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

std::string method_error(const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;
  os << std::format(R"(        throw std::runtime_error(
          std::format(R"-WSM-(Cannot execute method:

  {}()",
                    name);

  const auto               spaces      = std::string(name.size() + 1 + 2, ' ');
  const Size               largest_var = max_varlen(wsm);
  std::vector<std::string> arg;

  bool first = true;
  for (auto& t : wsm.out) {
    if (not first) os << ",\n" << spaces;
    first = false;
    std::print(os, R"({1} {0}: {{}})", std::string(largest_var - t.size(), ' '), t);
    arg.push_back(std::format(R"(_{0} ? "User-provided {1}"sv : "self.{0}"sv)", t, workspace_variables().at(t).type));
  }

  for (std::size_t i = 0; i < wsm.gout.size(); i++) {
    auto& t  = wsm.gout[i];
    auto& tt = wsm.gout_type[i];
    if (not first) os << ",\n" << spaces;
    first = false;
    std::print(os, R"({1} {0}: {{}})", std::string(largest_var - t.size(), ' '), t);
    if (tt == "Any") {
      arg.push_back(std::format(R"(_{0} ? std::format("User-provided {{}}", type(_{0})) : std::string("None"))", t));
    } else if (uses_variadic(tt)) {
      arg.push_back(std::format(R"(_{0} ? std::format("User-provided {{}}", type(_{0})) : std::string("None"))", t));
    } else {
      arg.push_back(std::format(R"(_{0} ? "User-provided {1}"sv : "self.{0}"sv)", t, tt));
    }
  }

  for (auto& t : wsm.in) {
    if (stdr::any_of(wsm.out, Cmp::eq(t))) continue;
    if (not first) os << ",\n" << spaces;
    first = false;
    std::print(os, R"({1} {0}: {{}})", std::string(largest_var - t.size(), ' '), t);
    arg.push_back(std::format(R"(_{0} ? "User-provided {1}"sv : "self.{0}"sv)", t, workspace_variables().at(t).type));
  }

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    auto& t  = wsm.gin[i];
    auto& v  = wsm.gin_value[i];
    auto& tt = wsm.gin_type[i];
    if (not first) os << ",\n" << spaces;
    first = false;
    std::print(os, R"({1} {0}: {{}})", std::string(largest_var - t.size(), ' '), t);

    if (tt == "Any") {
      arg.push_back(
          std::format(R"(_{0} ? std::format("User-provided {{}}", type(_{0})) : std::string(R"-WSMVAR-({1})-WSMVAR-"))",
                      t,
                      v ? std::format("{}", to_defval_str(*v, ""sv)) : "None"));
    } else if (uses_variadic(tt)) {
      arg.push_back(
          std::format(R"(_{0} ? std::format("User-provided {{}}", type(_{0})) : std::string(R"-WSMVAR-({1})-WSMVAR-"))",
                      t,
                      v ? std::format("{}", to_defval_str(*v, ""sv)) : "None"));
    } else {
      arg.push_back(
          std::format(R"(_{0} ? std::format("User-provided {{}}", type(_{0})) : std::string(R"-WSMVAR-({1})-WSMVAR-"))",
                      t,
                      v ? std::format("{}", to_defval_str(*v, ""sv)) : "None"));
    }
  }

  os << ")\n\nMethod reports the following error(s):\n{})-WSM-\", ";
  for (auto& a : arg) { os << std::format("{},\n           ", a); }
  os << "std::string_view(e.what())));\n";
  return os.str();
}

std::string method_argument_documentation(const WorkspaceMethodInternalRecord& wsm) {
  std::ostringstream os;

  bool first = true;

  for (const auto& t : wsm.out) {
    if (not first) os << ",\n    ";
    first = false;
    os << "\"" << t << "\"_a.noconvert().none() = py::none()";
  }

  for (Size i = 0; i < wsm.gout.size(); i++) {
    auto&  t  = wsm.gout[i];
    auto&& tt = wsm.gout_type[i];
    if (not first) os << ",\n    ";
    first = false;
    std::print(os, R"("{}"_a{}.none() = py::none())", t, tt == "Any" ? ""sv : ".noconvert()"sv);
  }

  for (const auto& t : wsm.in) {
    if (stdr::any_of(wsm.out, Cmp::eq(t))) continue;
    if (not first) os << ",\n    ";
    first = false;
    os << "\"" << t << "\"_a.none() = py::none()";
  }

  for (std::size_t i = 0; i < wsm.gin.size(); i++) {
    if (stdr::find(wsm.gout, wsm.gin[i]) != wsm.gout.end()) continue;
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

std::string method(const std::string& name, const WorkspaceMethodInternalRecord& wsm) {
  return std::format(
      R"-x-(  ws.def("{0}",[]({1}) -> {7} {{
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
      method_argument_selection(name, wsm) + method_input_checks(wsm),
      method_resolution(name, wsm),
      method_error(name, wsm),
      method_argument_documentation(wsm),
      method_docs(name),
      wsm.return_type);
}

void methods(int nfiles) {
  const auto& wsms = internal_workspace_methods();

  std::vector<std::ofstream> ofs(nfiles);
  for (int i = 0; i < nfiles; i++) { ofs[i].open(std::format("py_auto_wsm_{}.cpp", i)); }

  for (int i = 0; i < nfiles; i++) {
    select_ofstream(ofs, i) << R"--(#include <python_interface.h>

#include <workspace.h>

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

  for (int i = 0; i < nfiles; i++) { select_ofstream(ofs, i) << "}\n}  // namespace Python\n"; }
}
}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " output_file_count[int]" << '\n';
    return EXIT_FAILURE;
  }

  const int num_methods = std::stoi(argv[1]);

  methods(num_methods);

  if (errors.size()) {
    std::cerr << "Errors (" << errors.size() << "):\n";
    for (auto& e : errors) { std::cerr << e << '\n'; }
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
