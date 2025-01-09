#pragma once

#include <string>
#include <unordered_map>

struct WorkspaceVariableInternalRecord {
  std::string desc;
  std::string type;
  std::string default_value{""};
};

const std::unordered_map<std::string, WorkspaceVariableInternalRecord>&
internal_workspace_variables();

std::string_view any(const std::string& type);
