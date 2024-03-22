#pragma once

#include <string>
#include <unordered_map>

struct WorkspaceGroupRecord {
  std::string file;
  std::string desc;
  bool value_type{false};
};

const std::unordered_map<std::string, WorkspaceGroupRecord>&
internal_workspace_groups();
