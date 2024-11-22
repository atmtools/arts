#pragma once

#include <optional>
#include <string>
#include <unordered_map>

#include "workspace_wsv.h"

struct WorkspaceVariableInternalRecord {
  std::string desc;
  std::string type;
  std::optional<Wsv> default_value{std::nullopt};
};

std::unordered_map<std::string, WorkspaceVariableInternalRecord>
internal_workspace_variables();

std::string_view any(const std::string& type);
