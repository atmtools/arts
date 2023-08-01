#pragma once

#include <string>
#include <unordered_map>

struct WorkspaceGroupRecord {
  std::string file;
  std::string desc;
  WorkspaceGroupRecord() = default;
  WorkspaceGroupRecord(std::string&&, std::string&&);
};

std::unordered_map<std::string, WorkspaceGroupRecord> internal_workspace_groups();
