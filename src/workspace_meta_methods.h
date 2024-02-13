#pragma once

#include "workspace_methods.h"

struct WorkspaceMethodInternalMetaRecord {
  std::string name;
  std::string desc;
  std::vector<std::string> author;
  std::vector<std::string> methods;
  std::vector<std::string> out;

  [[nodiscard]] WorkspaceMethodInternalRecord create(
      const std::unordered_map<std::string, WorkspaceMethodInternalRecord>&
          wsms) const;
  [[nodiscard]] std::string call(
      const std::unordered_map<std::string, WorkspaceMethodInternalRecord>&
          wsms) const;
};

std::vector<WorkspaceMethodInternalMetaRecord> internal_meta_methods();
