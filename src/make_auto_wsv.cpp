#include <auto_wsg.h>
#include <debug.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "workspace_agendas.h"
#include "workspace_groups.h"
#include "workspace_variables.h"

const auto& wsg = internal_workspace_groups();
const auto& vars = internal_workspace_variables();
const auto& ags  = internal_workspace_agendas();

std::vector<std::string> wsv_names(std::ostream& os) {

  std::vector<std::string> names;
  names.reserve(vars.size() + ags.size());

  for (const auto& [name, data] : vars) {
    names.push_back(name);

    if (not wsg.contains(data.type)) {
      os << std::format(
          R"WSV(static_assert(false, R"--(

Workspace variable type is not a workspace group.

This is not allowed.

The workspace variable is:      "{}"
The workspace variable type is: "{}"
)--");

)WSV",
          name,
          data.type);
    }
  }

  for (const auto& [name, _] : ags) {
    names.push_back(name);
  }

  std::sort(names.begin(), names.end());

  auto ptr = std::adjacent_find(names.begin(), names.end());
  while (names.end() not_eq ptr) {
    os << std::format(R"WSV(static_assert(false, R"--(

Duplicate name in workspace variables and agendas.

This is not allowed.

The duplicate name is: "{}"
)--");

)WSV",
        *ptr);
    ptr = std::adjacent_find(ptr + 1, names.end());
  }

  return names;
}

std::vector<std::string> groups() {
  std::vector<std::string> groups;
  groups.reserve(wsg.size());

  for (const auto& [group_name, group] : wsg) {
    groups.push_back(group_name);
  }

  std::sort(groups.begin(), groups.end());

  return groups;
}

void header(std::ostream& os) {
  os << R"--(#pragma once

//! auto-generated by make_auto_wsv.cpp

#include <optional>
#include <string>
#include <unordered_map>

#include <auto_wsg.h>

#include <workspace_agenda_class.h>
#include <workspace_method_class.h>

struct WorkspaceVariableRecord {
  std::string desc;
  std::string type;
  std::optional<Wsv> default_value{};
};

const std::unordered_map<std::string, WorkspaceVariableRecord>& workspace_variables();
)--";
}

void implementation(std::ostream& os) {
  const static auto names = wsv_names(os);

  os << R"--(//! auto-generated by make_auto_wsv.cpp

#include "auto_wsv.h"

#include <workspace_agendas.h>
#include <workspace_variables.h>

std::unordered_map<std::string, WorkspaceVariableRecord> workspace_variables_create() {
  std::unordered_map<std::string, WorkspaceVariableRecord> vars;
  vars.reserve()--"
     << names.size() << R"--();

  for (auto&& [name, record] : internal_workspace_variables()) {
    vars.emplace(name, WorkspaceVariableRecord{record.desc, record.type, record.default_value});
  }

  for (auto&& [name, record] : internal_workspace_agendas()) {
    vars.emplace(name, WorkspaceVariableRecord{record.desc, record.array ? "ArrayOfAgenda" : "Agenda"});
  }

  return vars;
}

const std::unordered_map<std::string, WorkspaceVariableRecord>& workspace_variables() {
  const static auto vars = workspace_variables_create();
  return vars;
}
)--";
}

int main() try {
  std::ofstream head("auto_wsv.h");
  std::ofstream impl("auto_wsv.cpp");

  header(head);
  implementation(impl);
} catch (std::exception& e) {
  std::cerr << "Cannot create the automatic variables with error:\n\n"
            << e.what() << '\n';
  return 1;
}
