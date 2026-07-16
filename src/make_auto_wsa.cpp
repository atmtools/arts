#include <array_algo.h>
#include <auto_wsg.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <map>
#include <ranges>
#include <utility>

#include "workspace_agendas.h"
#include "workspace_variables.h"

namespace {
const auto& wsa = internal_workspace_agendas();
const auto& wsv = internal_workspace_variables();

//! List of groups and their variables
struct auto_ag {
  std::string                                      desc;
  std::vector<std::pair<std::string, std::string>> o;
  std::vector<std::pair<std::string, std::string>> i;
  std::vector<StringVectorAgendaHelper>            output_constraints;
};

void helper_auto_ag(std::ostream&                                     os,
                    std::vector<std::pair<std::string, std::string>>& vars,
                    const std::string&                                name,
                    const std::string&                                var) {
  auto ptr = wsv.find(var);
  if (ptr == wsv.end()) {
    auto ag_ptr = wsa.find(var);
    if (ag_ptr == wsa.end()) {
      std::print(os,
                 R"X(static_assert(false, R"-err-(

Could not find workspace variable "{}" of agenda "{}"
)-err-");)X",
                 var,
                 name);
    } else {
      vars.emplace_back("Agenda", var);
    }
  } else {
    vars.emplace_back(ptr->second.type, var);
  }
}

std::map<std::string, auto_ag> auto_ags(std::ostream& os) {
  std::map<std::string, auto_ag> map;

  for (const auto& [name, record] : wsa) {
    auto& ag = map[name];
    ag.desc  = record.desc;

    for (const auto& out : record.output) { helper_auto_ag(os, ag.o, name, out); }

    for (const auto& in : record.input) { helper_auto_ag(os, ag.i, name, in); }

    ag.output_constraints = record.output_constraints;
  }

  return map;
}

void header_docstring(std::ostream& os, const auto_ag& ag, const std::string& agname) {
  std::println(os, "/** {}\n  @param[in] ws The workspace", ag.desc);

  const auto is_input = [&ag](auto& v) {
    return std::find_if(ag.i.begin(), ag.i.end(), [v](auto& i) { return i.second == v; }) != ag.i.end();
  };
  const auto is_output = [&ag](auto& v) {
    return std::find_if(ag.o.begin(), ag.o.end(), [v](auto& o) { return o.second == v; }) != ag.o.end();
  };

  for (auto& [type, name] : ag.o) {
    std::println(os, "  @param[{}] {} As WSV", is_input(name) ? "inout"sv : "out"sv, name);
  }

  for (auto& [type, name] : ag.i) {
    if (not is_output(name)) std::println(os, "  @param[in] {} As WSV", name);
  }

  std::println(os,
               R"(  @param[in] {} As WSV
  @param[in] local_workspace_pointer - will be moved from if given
  @return Workspace - the shared and copied values of all the work the agenda performed
*/)",
               agname);
}

void call_operator(std::ostream& os, const auto_ag& ag, const std::string& agname, bool is_header) {
  const auto is_output = [&ag](auto& v) {
    return std::find_if(ag.o.begin(), ag.o.end(), [v](auto& o) { return o.second == v; }) != ag.o.end();
  };

  const std::string spaces(18 + agname.size(), ' ');

  std::print(os, "Workspace {}Execute(const Workspace& ws", agname);

  for (auto& [type, name] : ag.o) { std::print(os, ",\n{}{}& {}", spaces, type, name); }

  for (auto& [type, name] : ag.i) {
    if (not is_output(name)) std::print(os, ",\n{}const {}& {}", spaces, type, name);
  }

  std::print(
      os, ",\n{0}const Agenda& {1},\n{0}Workspace* _lws_ptr{2})", spaces, agname, is_header ? " = nullptr"sv : ""sv);
}

void header(std::ostream& os) {
  const auto agmap = auto_ags(os);

  std::println(os, R"--(#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "auto_wsg.h"

struct Workspace;

struct WorkspaceAgendaRecord {{
  std::string desc;
  std::vector <std::string> output;
  std::vector <std::string> input;
  }};

const std::unordered_map<std::string, WorkspaceAgendaRecord>& workspace_agendas();

// Returns documentation of default enums as keys for agenda names and then a list of pairs of option and documentation
[[nodiscard]]
std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> get_agenda_enum_documentation();

struct WorkspaceAgendaBoolHandler {{
)--");

  for (auto& ag : agmap) { std::println(os, "  bool has_{} : 1 {{false}};", ag.first); }

  std::println(os, R"--(
  [[nodiscard]] bool has(const std::string&) const;
  void set(const std::string&);
  friend std::ostream& operator<<(std::ostream& os, const WorkspaceAgendaBoolHandler& wab);
}};
)--");

  for (auto& ag : agmap) {
    header_docstring(os, ag.second, ag.first);
    call_operator(os, ag.second, ag.first, true);
    std::println(os, ";\n");
  }

  for (auto& [name, ag] : wsa) {
    if (ag.enum_options.empty()) continue;

    std::println(os,
                 R"(
void {0}Set(Agenda& {0}, const String&);
Agenda get_{0}(const std::string_view);
)",
                 name);
  }
}

void agenda_checker(std::ostream& os, const std::string& name) {
  std::println(os,
               R"xx(  if (not {0}.is_checked()) {{
    throw std::runtime_error(R"--(
You have somehow created the agenda without checking it.

Please manually call finalize() on the agenda
)--");
  }}
)xx",
               name);
}

std::string double_curly(std::string s) {
  Size n = 0;

  for (Size i = 0; i < s.size() - 1; i++) { n += (s[i] == '{' and s[i + 1] == '{'); }

  for (Size i = 0; i < n; i++) {
    s.push_back('{');
    s.push_back('}');
  }

  for (Size i = 1; i < s.size(); i++) {
    if (s[i - 1] == '{' and s[i] == '}') { stdr::rotate(s.begin() + i, s.end() - 2, s.end()); }
  }

  return s;
}

void workspace_setup_and_exec(std::ostream& os, const std::string& name, const auto_ag& ag) {
  std::println(os, R"(
  // Create a local workspace upon need or get the data from the pointer
  Workspace _lws = _lws_ptr ? std::move(*_lws_ptr) : Workspace(WorkspaceInitialization::Empty);
  const bool empty = _lws.wsv.empty();

  // Name the original data here)");
  for (auto& i : ag.i) { std::println(os, R"(  static const std::string _wsv_{0} = "{0}";)", i.second); }

  for (auto& o : ag.o) {
    if (stdr::find_if(ag.i, [&o](auto& v) { return v.second == o.second; }) != ag.i.end()) continue;
    std::println(os, R"(  static const std::string _wsv_{0} = "{0}";)", o.second);
  }

  std::println(os, R"(
  // Always share original data here)");
  for (auto& i : ag.i) { std::println(os, R"(  _lws.set(_wsv_{1}, const_cast<{0}*>(&{1}));)", i.first, i.second); }

  std::println(os,
               R"(
  // Copy and share data from old workspace (this will copy pure inputs that are modified)
  if (empty)
    {0}.copy_workspace(_lws, ws);
  else
    {0}.copy_only_workspace(_lws, ws);

  // Modified data must be copied here)",
               name);
  for (auto& o : ag.o) {
    std::println(os, R"(  _lws.overwrite(_wsv_{1}, const_cast<{0}*>(&{1}));)", o.first, o.second);
  }

  std::println(os, "\n  // Run all the methods\n  {}.execute(_lws);", name);

  for (auto& constraint : ag.output_constraints) {
    std::println(os,
                 R"--(
  if(not ({}))
    throw std::runtime_error(std::format(R"ERR({})--",
                 constraint.test,
                 double_curly(constraint.constraint));

    const Size N = max(constraint.printables, [](const std::string& x) -> Size { return x.size(); });
    for (auto& p : constraint.printables) {
      std::println(os,
                   R"--(
{}:{} {{}})--",
                   double_curly(p),
                   std::string(N - p.size(), ' '));
    }
    if (not constraint.printables.empty()) {
      std::print(os,
                 R"--(
)ERR", {:,}));)--",
                 constraint.printables);
    }
  }

  std::println(os, R"(
  // Remove the unsafe content (false sharing pointers))");
  for (auto& i : ag.i) { std::println(os, R"(  _lws.erase(_wsv_{1});)", i.first, i.second); }

  for (auto& i : ag.o) {
    if (stdr::find_if(ag.i, [&i](auto& v) { return v.second == i.second; }) != ag.i.end()) continue;
    std::println(os, R"(  _lws.erase(_wsv_{1});)", i.first, i.second);
  }

  std::println(os, "\n  return _lws;");
}

void implementation(std::ostream& os) {
  const auto agmap = auto_ags(os);

  std::println(os,
               R"--(#include "auto_wsa.h"

#include <workspace_agenda_class.h>
#include <workspace_method_class.h>
#include <workspace_agenda_creator.h>
#include <workspace_class.h>
#include <time_report.h>

std::unordered_map<std::string, WorkspaceAgendaRecord> get_workspace_agendas() {{
  std::unordered_map<std::string, WorkspaceAgendaRecord> ags;
  ags.reserve({});
)--",
               wsa.size());

  auto quoted = std::views::transform([](const std::string& s) { return '"' + s + '"'; });

  for (const auto& [name, ag] : agmap) {
    std::println(os,
                 R"-x-(ags["{}"] = WorkspaceAgendaRecord{{
    .desc=R"--({})--",
    .output={{{:,}}},
    .input={{{:,}}},
}};
)-x-",
                 name,
                 ag.desc,
                 ag.o | stdv::values | quoted | stdr::to<std::vector<std::string>>(),
                 ag.i | stdv::values | quoted | stdr::to<std::vector<std::string>>());
  }

  std::println(os, R"--(  return ags;
}}

const std::unordered_map<std::string, WorkspaceAgendaRecord>& workspace_agendas() {{
  const static auto ags = get_workspace_agendas();
  return ags;
}}

bool WorkspaceAgendaBoolHandler::has(const std::string& ag) const {{)--");

  for (auto& ag : agmap) { std::println(os, R"(  if (ag == "{0}") return has_{0};)", ag.first); }
  std::println(os, R"--(
  throw std::runtime_error(std::format("Not a predefined agenda: \"{{}}\"", ag));
}}

void WorkspaceAgendaBoolHandler::set(const std::string& ag) {{
)--");
  for (auto& ag : agmap) { std::println(os, R"(  if (ag == "{0}") {{has_{0} = true; return;}})", ag.first, ag.first); }
  std::println(os, R"--(
  throw std::runtime_error(std::format("Not a predefined agenda: \"{{}}\"", ag));
}}
std::ostream& operator<<(std::ostream& os, const WorkspaceAgendaBoolHandler& wab) {{
)--");

  const Size n =
      stdr::max_element(agmap | stdv::keys, {}, [](const std::string& s) { return s.size(); }).base()->first.size();
  for (auto& ag : agmap) {
    std::println(
        os, R"(  os << "{0}: {2}" << wab.has_{0} << '\n';)", ag.first, ag.first, std::string(n - ag.first.size(), ' '));
  }
  std::println(os, R"--(
  return os;
}}
)--");
  for (const auto& [name, ag] : agmap) {
    call_operator(os, ag, name, false);
    std::println(os, " try {{\n  ARTS_TIME_REPORT\n");
    agenda_checker(os, name);
    workspace_setup_and_exec(os, name, ag);
    std::println(os,
                 R"(}} catch(std::exception& e) {{
  throw std::runtime_error(std::format(R"--(Error executing agenda "{}":
{{}})--", e.what()));
}}
)",
                 name);
  }

  for (auto& [name, ag] : wsa) {
    if (ag.enum_options.empty()) continue;

    std::print(os,
               R"(
void {0}Set(Agenda& ag, const String& option) try {{
  ARTS_TIME_REPORT

  ag = get_{0}(option);  
}}  ARTS_METHOD_ERROR_CATCH
)",
               name);
  }

  std::print(os, R"(
std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> get_agenda_enum_documentation() {{
  std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> out{{}};

)");

  for (auto& [name, ag] : wsa) {
    std::print(os,
               R"(  [[maybe_unused]] auto& tmp_{0} = out["{0}Set"];

)",
               name);
    if (ag.enum_options.empty()) continue;

    for (auto& opt : ag.enum_options) {
      std::print(os,
                 R"(  tmp_{1}.emplace_back("{0}", get_{1}("{0}").sphinx_list());
)",
                 opt,
                 name);
    }
  }

  std::println(os, "  return out;\n}}");
}

void options(std::ostream& os) {
  std::print(os, R"(#pragma once

#include <enums.h>

)");

  for (auto& [name, ag] : wsa) {
    if (ag.enum_options.empty()) continue;

    std::print(os, "enum class {}Predefined {{\n", name);
    for (auto& opt : ag.enum_options) { std::print(os, "  {},\n", opt); }

    std::print(os,
               R"(}};

template <>
constexpr {0}Predefined to<{0}Predefined>(const std::string_view x) {{
)",
               name);

    for (auto& opt : ag.enum_options) {
      std::print(os,
                 R"(  if (x == "{0}"sv) return {1}Predefined::{0};
)",
                 opt,
                 name);
    }

    std::print(os,
               R"-XX-(  throw std::runtime_error(R"-X-(Bad value. Valid options for "{0}":

{1:BNq,}
)-X-");
}}

)-XX-",
               name,
               ag.enum_options);
  }
}
}  // namespace

int main() try {
  std::ofstream head("auto_wsa.h");
  std::ofstream impl("auto_wsa.cpp");
  std::ofstream optshh("auto_wsa_options.h");

  header(head);
  implementation(impl);
  options(optshh);
} catch (std::exception& e) {
  std::println(stderr, "Cannot create the automatic agendas with error:\n\n{}", e.what());
  return 1;
}
