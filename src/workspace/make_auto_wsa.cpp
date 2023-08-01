#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <ostream>
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

void helper_auto_ag(std::vector<std::pair<std::string, std::string>>& vars, const std::string& name, const std::string& var) {
  auto ptr = wsv.find(var);
  if (ptr == wsv.end()) {
    auto ag_ptr = wsa.find(var);
    if (ag_ptr == wsa.end()) {
      std::cerr << R"(static_assert(false, "\n\nCould not find workspace variable \")" << var << R"(\" of agenda \")" << name << "\\\"\\n\\n\");\n";
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

std::map<std::string, auto_ag> auto_ags() {
  std::map<std::string, auto_ag> map;

  for (const auto& [name, record] : wsa) {
    auto& ag = map[name];
    ag.desc = record.desc;
    ag.array = record.array;

    for (const auto& out: record.output) {
      helper_auto_ag(ag.o, name, out);
    }

    for (const auto& in: record.input) {
      if (std::none_of(record.output.begin(), record.output.end(), [in](auto& elem){return elem == in;})) {
        helper_auto_ag(ag.i, name, in);
      }
    }
  }

  return map;
}

void header_docstring(std::ostream& os,
                      const auto_ag& ag,
                      const std::string& agname) {
  os << "/** " << ag.desc << "\n  @param[in] ws The workspace";

  for (auto& [type, name] : ag.o) {
    if (std::find_if(ag.i.begin(), ag.i.end(), [name](auto& v){return v.second == name;}) == ag.i.end())
      os << "\n  @param[out] " << name << " As WSV";
    else
      os << "\n  @param[inout] " << name << " As WSV";
  }

  for (auto& [type, name] : ag.i) {
    os << "\n  @param[in] " << name << " As WSV";
  }

  os << "\n  @param[in] " << agname << " As WSV"  << "\n*/\n";
}

void call_operator(std::ostream& os,
                   const auto_ag& ag,
                   const std::string& agname) {
  const std::string spaces(5 + agname.size() + 8, ' ');

  os << "void " << agname << "Execute(const Workspace& ws";

  for (auto& [type, name] : ag.o) {
    os << ",\n" << spaces << type << "& " << name;
  }

  for (auto& [type, name] : ag.i) {
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
  const auto agmap = auto_ags();

  os << R"--(#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include <auto_wsg.h>

struct WorkspaceAgendaRecord {
  std::string desc;
  std::vector <std::string> inout;
  std::vector <std::string> output;
  std::vector <std::string> input;
};

const std::unordered_map<std::string, WorkspaceAgendaRecord>& workspace_agendas();

)--";

  for (auto& ag: agmap) {
    header_docstring(os, ag.second, ag.first);
    call_operator(os, ag.second, ag.first);
    os << ";\n\n";
  }
}

void agenda_checker(std::ostream& os, const std::string& name, bool array) {
  if (array) {
    os << "  if (agenda_array_index < 0 or static_cast<std::size_t>(agenda_array_index) >= " << name << ".size()) {\n"
          "    throw std::runtime_error(R\"--(Agenda "
       << std::quoted(name) << " has invalid index)--\");\n  }\n\n";
  }

  os << "  if (not " << name;
  if (array) {
    os << "[agenda_array_index]";
  }
  os << ".is_checked()) {\n"
        "    throw std::runtime_error(R\"--(Agenda "
     << std::quoted(name) << " is not checked)--\");\n  }\n";
}

void workspace_setup_and_exec(std::ostream& os, const std::string& name, const auto_ag& ag) {
  os << "\n  Workspace _lws{WorkspaceInitialization::Empty};\n";

  for (auto& o: ag.o) {
    os << "  _lws.set(" << std::quoted(o.second) << ", &" << o.second << ");\n";
  }

  for (auto& i: ag.i) {
    os << "  _lws.set(" << std::quoted(i.second) << ", std::make_shared<Wsv>(" << i.second << "));\n";
  }

  os << "  " << name;
  if (ag.array) {
    os << "[agenda_array_index]";
  }
  os << ".copy_workspace(_lws, ws);\n\n";

  os << "  " << name;
  if (ag.array) {
    os << "[agenda_array_index]";
  }
  os << ".execute(_lws);\n";
}

void implementation(std::ostream& os) {
  const auto agmap = auto_ags();

  os << R"--(#include <auto_wsa.h>

#include "workspace_agenda_class.h"
#include "workspace_method_class.h"

#include "workspace_class.h"

std::unordered_map<std::string, WorkspaceAgendaRecord> get_workspace_agendas() {
  std::unordered_map<std::string, WorkspaceAgendaRecord> ags;
  ags.reserve()--"
     << wsa.size() << ");\n\n";

  for (const auto& [name, ag] : agmap) {
    os << "ags[" << std::quoted(name)
       << "] = WorkspaceAgendaRecord{\n    R\"--(" << ag.desc << ")--\","
       << "\n    {";
    os << "},\n    {";
    for (const auto& o : ag.o) {
      os << std::quoted(o.second) << ", ";
    }
    os << "},\n    {";
    for (const auto& i : ag.i) {
      os << std::quoted(i.second) << ", ";
    }
    os << "}\n  };\n\n";
  }

  os << R"--(  return ags;
}

const std::unordered_map<std::string, WorkspaceAgendaRecord>& workspace_agenda() {
  const static auto ags = get_workspace_agendas();
  return ags;
}

)--";
  for (const auto& [name, ag] : agmap) {
    call_operator(os, ag, name);
    os << " try {\n";
    agenda_checker(os, name, ag.array);
    workspace_setup_and_exec(os, name, ag);
    os << "} catch(std::exception& e) {\n  throw std::runtime_error(var_string(R\"--(Error executing agenda "
       << std::quoted(name) << ":\n)--\", e.what()));\n}\n\n";
  }
}

int main() {
header(std::cout);
implementation(std::cerr);
}
