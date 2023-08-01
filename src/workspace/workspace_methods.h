#pragma once

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <auto_wsg.h>

struct WorkspaceMethodInternalRecord {
  std::string  desc{};
  std::vector<std::string> author{};
  std::vector<std::string> out{};
  std::vector<std::string> gout{};
  std::vector<std::string> gout_type{};
  std::vector<std::string> gout_desc{};
  std::vector<std::string> in{};
  std::vector<std::string> gin{};
  std::vector<std::string> gin_type{};
  std::vector<std::optional<Wsv>> gin_value{};
  std::vector<std::string> gin_desc{};
  bool pass_workspace{false};
  bool pass_names{false};
};

std::unordered_map<std::string, WorkspaceMethodInternalRecord> internal_workspace_methods();
