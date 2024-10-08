#include <algorithm>
#include <exception>
#include <iomanip>
#include <iostream>
#include <map>
#include <utility>

#include "workspace_agendas.h"
#include "workspace_variables.h"

const auto wsa = internal_workspace_agendas();
const auto wsv = internal_workspace_variables();

//! List of groups and their variables
struct auto_ag {
  std::string desc;
  std::vector<std::pair<std::string, std::string>> o;
  std::vector<std::pair<std::string, std::string>> i;
  bool array;
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
      if (ag_ptr->second.array) {
        vars.emplace_back("ArrayOfAgenda", var);
      } else {
        vars.emplace_back("Agenda", var);
      }
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
    ag.array = record.array;

    for (const auto& out : record.output) {
      helper_auto_ag(os, ag.o, name, out);
    }

    for (const auto& in : record.input) {
      helper_auto_ag(os, ag.i, name, in);
    }
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

  os << ",\n" << spaces << "const ";

  if (ag.array) {
    os << "ArrayOf";
  }

  os << "Agenda& " << agname;
  os << ")";
}

void header(std::ostream& os) {
  const auto agmap = auto_ags(os);

  os << R"--(#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include <auto_wsg.h>

class Workspace;

struct WorkspaceAgendaRecord {
  std::string desc;
  std::vector <std::string> output;
  std::vector <std::string> input;
};

const std::unordered_map<std::string, WorkspaceAgendaRecord>& workspace_agendas();

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
}

void agenda_checker(std::ostream& os, const std::string& name, bool array) {
  if (array) {
    os << "  if (agenda_array_index < 0 or static_cast<std::size_t>(agenda_array_index) >= "
       << name
       << ".size()) {\n"
          "    throw std::runtime_error(R\"--(Array index out-of-bounds)--\");\n  }\n\n";
  }

  os << "  if (not " << name;
  if (array) {
    os << "[agenda_array_index]";
  }
  os << ".is_checked()) {\n"
        "    throw std::runtime_error(R\"--(\nYou have somehow created the agenda without checking it.\n\nPlease manually call finalize() on the agenda)--\");\n  }\n\n";

  os << "  if (const auto& n = " << name;
  if (array) {
    os << "[agenda_array_index]";
  }
  os << ".get_name(); n != \"" << name
     << "\") {\n"
        "    throw std::runtime_error(var_string(\"Mismatch with name: \", n));\n  }\n";
}

void workspace_setup_and_exec(std::ostream& os,
                              const std::string& name,
                              const auto_ag& ag) {
  os << "\n  Workspace _lws{WorkspaceInitialization::Empty};\n\n";

  os << "  // Always share original data here\n";

  // FIXME: This should be overwrite, no?  And then if it exists we should copy over it?
  for (auto& i : ag.i) {
    os << "  _lws.set(" << var_string('"', i.second, '"') << ", const_cast<"
       << i.first << "*>(&" << i.second << "));\n";
  }

  os << "\n"
        "  // Copy and share data from old workspace (this will copy pure inputs that are modified)\n  ";
  os << name;
  if (ag.array) {
    os << "[agenda_array_index]";
  }
  os << ".copy_workspace(_lws, ws);\n";

  os << "\n  // Modified data must be copied here\n";
  for (auto& o : ag.o) {
    os << "  _lws.overwrite(" << var_string('"', o.second, '"') << ", &"
       << o.second << ");\n";
  }

  os << "\n  // Run all the methods\n  " << name;
  if (ag.array) {
    os << "[agenda_array_index]";
  }
  os << ".execute(_lws);\n";
}

void implementation(std::ostream& os) {
  const auto agmap = auto_ags(os);

  os << R"--(#include <auto_wsa.h>

#include "workspace_agenda_class.h"
#include "workspace_method_class.h"

#include "workspace_class.h"

std::unordered_map<std::string, WorkspaceAgendaRecord> get_workspace_agendas() {
  std::unordered_map<std::string, WorkspaceAgendaRecord> ags;
  ags.reserve()--"
     << wsa.size() << ");\n\n";

  for (const auto& [name, ag] : agmap) {
    os << "ags[\"" << name << '"'
       << "] = WorkspaceAgendaRecord{\n    .desc=R\"--(" << ag.desc << ")--\","
       << "\n    .output={";
    for (const auto& o : ag.o) {
      os << var_string('"', o.second, '"') << ", ";
    }
    os << "},\n    .input={";
    for (const auto& i : ag.i) {
      os << var_string('"', i.second, '"') << ", ";
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
  os << R"--(
  throw std::runtime_error(var_string("Not a predefined agenda: \"", ag, '"'));
}

void WorkspaceAgendaBoolHandler::set(const std::string& ag) {
)--";
  for (auto& ag : agmap) {
    os << "  if (ag == \"" << ag.first << "\") {has_" << ag.first
       << " = true; return;}\n";
  }
  os << R"--(
  throw std::runtime_error(var_string("Not a predefined agenda: \"", ag, '"'));
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
    os << " try {\n";
    agenda_checker(os, name, ag.array);
    workspace_setup_and_exec(os, name, ag);
    os << "} catch(std::exception& e) {\n  throw std::runtime_error(std::format(R\"--(Error executing agenda "
       << '"' << name << '"' << ":\n{})--\", e.what()));\n}\n\n";
  }
}

int main() try {
  std::ofstream head("auto_wsa.h");
  std::ofstream impl("auto_wsa.cpp");

  header(head);
  implementation(impl);
} catch (std::exception& e) {
  std::cerr << "Cannot create the automatic agendas with error:\n\n"
            << e.what() << '\n';
  return 1;
}
