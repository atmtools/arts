#include <array_algo.h>
#include <auto_wsg.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <map>
#include <utility>

#include "workspace_agendas.h"
#include "workspace_variables.h"

const auto& wsa = internal_workspace_agendas();
const auto& wsv = internal_workspace_variables();

//! List of groups and their variables
struct auto_ag {
  std::string desc;
  std::vector<std::pair<std::string, std::string>> o;
  std::vector<std::pair<std::string, std::string>> i;
  std::vector<StringVectorAgendaHelper> output_constraints;
};

void helper_auto_ag(std::ostream& os,
                    std::vector<std::pair<std::string, std::string>>& vars,
                    const std::string& name,
                    const std::string& var) {
  auto ptr = wsv.find(var);
  if (ptr == wsv.end()) {
    auto ag_ptr = wsa.find(var);
    if (ag_ptr == wsa.end()) {
      os << R"(static_assert(false, "\n\nCould not find workspace variable \")"
         << var << R"(\" of agenda \")" << name << "\\\"\\n\\n\");\n";
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

    for (const auto& out : record.output) {
      helper_auto_ag(os, ag.o, name, out);
    }

    for (const auto& in : record.input) {
      helper_auto_ag(os, ag.i, name, in);
    }

    ag.output_constraints = record.output_constraints;
  }

  return map;
}

void header_docstring(std::ostream& os,
                      const auto_ag& ag,
                      const std::string& agname) {
  os << "/** " << ag.desc << "\n  @param[in] ws The workspace";

  const auto is_input = [&ag](auto& v) {
    return std::find_if(ag.i.begin(), ag.i.end(), [v](auto& i) {
             return i.second == v;
           }) != ag.i.end();
  };
  const auto is_output = [&ag](auto& v) {
    return std::find_if(ag.o.begin(), ag.o.end(), [v](auto& o) {
             return o.second == v;
           }) != ag.o.end();
  };

  for (auto& [type, name] : ag.o) {
    if (not is_input(name))
      os << "\n  @param[out] " << name << " As WSV";
    else
      os << "\n  @param[inout] " << name << " As WSV";
  }

  for (auto& [type, name] : ag.i) {
    if (not is_output(name)) os << "\n  @param[in] " << name << " As WSV";
  }

  os << "\n  @param[in] " << agname << " As WSV"
     << "\n*/\n";
}

void call_operator(std::ostream& os,
                   const auto_ag& ag,
                   const std::string& agname) {
  const auto is_output = [&ag](auto& v) {
    return std::find_if(ag.o.begin(), ag.o.end(), [v](auto& o) {
             return o.second == v;
           }) != ag.o.end();
  };

  const std::string spaces(5 + agname.size() + 8, ' ');

  os << "void " << agname << "Execute(const Workspace& ws";

  for (auto& [type, name] : ag.o) {
    os << ",\n" << spaces << type << "& " << name;
  }

  for (auto& [type, name] : ag.i) {
    if (not is_output(name))
      os << ",\n" << spaces << "const " << type << "& " << name;
  }

  os << ",\n" << spaces << "const Agenda& " << agname << ")";
}

void header(std::ostream& os) {
  const auto agmap = auto_ags(os);

  os << R"--(#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "auto_wsg.h"

struct Workspace;

struct WorkspaceAgendaRecord {
  std::string desc;
  std::vector <std::string> output;
  std::vector <std::string> input;
};

const std::unordered_map<std::string, WorkspaceAgendaRecord>& workspace_agendas();

// Returns documentation of default enums as keys for agenda names and then a list of pairs of option and documentation
[[nodiscard]]
std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> get_agenda_enum_documentation();

struct WorkspaceAgendaBoolHandler {
)--";

  for (auto& ag : agmap) {
    os << "  bool has_" << ag.first << " : 1 {false};\n";
  }

  os << R"--(
  [[nodiscard]] bool has(const std::string&) const;
  void set(const std::string&);
  friend std::ostream& operator<<(std::ostream& os, const WorkspaceAgendaBoolHandler& wab);
};

)--";

  for (auto& ag : agmap) {
    header_docstring(os, ag.second, ag.first);
    call_operator(os, ag.second, ag.first);
    os << ";\n\n";
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
  os << "  if (not " << name;
  os << ".is_checked()) {\n"
        "    throw std::runtime_error(R\"--(\nYou have somehow created the agenda without checking it.\n\nPlease manually call finalize() on the agenda)--\");\n  }\n\n";

  os << "  if (const auto& n = " << name;
  os << ".get_name(); n != \"" << name
     << "\") {\n"
        "    throw std::runtime_error(std::format(\"Mismatch with name: {}\", n));\n  }\n";
}

std::string double_curly(std::string s) {
  Size n = 0;

  for (Size i = 0; i < s.size() - 1; i++) {
    n += (s[i] == '{' and s[i + 1] == '{');
  }

  for (Size i = 0; i < n; i++) {
    s.push_back('{');
    s.push_back('}');
  }

  for (Size i = 1; i < s.size(); i++) {
    if (s[i - 1] == '{' and s[i] == '}') {
      stdr::rotate(s.begin() + i, s.end() - 2, s.end());
    }
  }

  return s;
}

void workspace_setup_and_exec(std::ostream& os,
                              const std::string& name,
                              const auto_ag& ag) {
  os << "\n  Workspace _lws{WorkspaceInitialization::Empty};\n\n";

  os << "  // Always share original data here\n";

  // FIXME: This should be overwrite, no?  And then if it exists we should copy over it?
  for (auto& i : ag.i) {
    os << "  _lws.set(" << std::format(R"("{}")", i.second) << ", const_cast<"
       << i.first << "*>(&" << i.second << "));\n";
  }

  os << "\n"
        "  // Copy and share data from old workspace (this will copy pure inputs that are modified)\n  ";
  os << name;
  os << ".copy_workspace(_lws, ws);\n";

  os << "\n  // Modified data must be copied here\n";
  for (auto& o : ag.o) {
    os << "  _lws.overwrite(" << std::format(R"("{}")", o.second) << ", &"
       << o.second << ");\n";
  }

  os << "\n  // Run all the methods\n  " << name;
  os << ".execute(_lws);\n";

  for (auto& constraint : ag.output_constraints) {
    std::print(os,
               R"--(
  if(not ({}))
    throw std::runtime_error(std::format(R"ERR({}
)--",
               constraint.test,
               double_curly(constraint.constraint));

    const Size N = max(constraint.printables,
                       [](const std::string& x) -> Size { return x.size(); });
    for (auto& p : constraint.printables) {
      std::print(os,
                 R"--(
{}:{} {{}})--",
                 double_curly(p),
                 std::string(N - p.size(), ' '));
    }
    if (not constraint.printables.empty()) {
      std::print(os,
                 R"--(
)ERR", {:,}));
)--",
                 constraint.printables);
    }
  }
}

void implementation(std::ostream& os) {
  const auto agmap = auto_ags(os);

  os << R"--(#include "auto_wsa.h"

#include <workspace_agenda_class.h>
#include <workspace_method_class.h>
#include <workspace_agenda_creator.h>
#include <workspace_class.h>
#include <time_report.h>

std::unordered_map<std::string, WorkspaceAgendaRecord> get_workspace_agendas() {
  std::unordered_map<std::string, WorkspaceAgendaRecord> ags;
  ags.reserve()--"
     << wsa.size() << ");\n\n";

  for (const auto& [name, ag] : agmap) {
    os << "ags[\"" << name << '"'
       << "] = WorkspaceAgendaRecord{\n    .desc=R\"--(" << ag.desc << ")--\","
       << "\n    .output={";
    for (const auto& o : ag.o) {
      os << std::format(R"("{}")", o.second) << ", ";
    }
    os << "},\n    .input={";
    for (const auto& i : ag.i) {
      os << std::format(R"("{}")", i.second) << ", ";
    }
    os << "}\n  };\n\n";
  }

  os << R"--(  return ags;
}

const std::unordered_map<std::string, WorkspaceAgendaRecord>& workspace_agendas() {
  const static auto ags = get_workspace_agendas();
  return ags;
}

bool WorkspaceAgendaBoolHandler::has(const std::string& ag) const {
)--";
  for (auto& ag : agmap) {
    os << "  if (ag == \"" << ag.first << "\") return has_" << ag.first
       << ";\n";
  }
  for (auto& ag : internal_workspace_agenda_names()) {
    os << "  if (ag == \"" << ag.first << "\") return has_" << ag.second
       << ";\n";
  }
  os << R"--(
  throw std::runtime_error(std::format("Not a predefined agenda: \"{}\"", ag));
}

void WorkspaceAgendaBoolHandler::set(const std::string& ag) {
)--";
  for (auto& ag : agmap) {
    os << "  if (ag == \"" << ag.first << "\") {has_" << ag.first
       << " = true; return;}\n";
  }
  for (auto& ag : internal_workspace_agenda_names()) {
    os << "  if (ag == \"" << ag.first << "\") {has_" << ag.second
       << " = true; return;}\n";
  }
  os << R"--(
  throw std::runtime_error(std::format("Not a predefined agenda: \"{}\"", ag));
}
std::ostream& operator<<(std::ostream& os, const WorkspaceAgendaBoolHandler& wab) {
)--";
  for (auto& ag : agmap) {
    os << "  os << \"" << ag.first << ": \" << wab.has_" << ag.first
       << " << '\\n';\n";
  }
  os << R"--(
  return os;
}

)--";
  for (const auto& [name, ag] : agmap) {
    call_operator(os, ag, name);
    os << " try {\n  ARTS_TIME_REPORT\n\n";
    agenda_checker(os, name);
    workspace_setup_and_exec(os, name, ag);
    os << "} catch(std::exception& e) {\n  throw std::runtime_error(std::format(R\"--(Error executing agenda "
       << '"' << name << '"' << ":\n{})--\", e.what()));\n}\n\n";
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

    os << '\n';
  }

  os << "  return out;\n}\n";
}

void options(std::ostream& os) {
  std::print(os, R"(#pragma once

#include <enums.h>

)");

  for (auto& [name, ag] : wsa) {
    if (ag.enum_options.empty()) continue;

    std::print(os, "enum class {}Predefined {{\n", name);
    for (auto& opt : ag.enum_options) {
      std::print(os, "  {},\n", opt);
    }

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

    std::print(
        os,
        R"-XX-(  throw std::runtime_error(R"-X-(Bad value. Valid options for "{0}":

{1:BNq,}
)-X-");
}}

)-XX-",
        name,
        ag.enum_options);
  }
}

int main() try {
  std::ofstream head("auto_wsa.h");
  std::ofstream impl("auto_wsa.cpp");
  std::ofstream optshh("auto_wsa_options.h");

  header(head);
  implementation(impl);
  options(optshh);
} catch (std::exception& e) {
  std::cerr << "Cannot create the automatic agendas with error:\n\n"
            << e.what() << '\n';
  return 1;
}
