#pragma once

#include <auto_wsv_value_wrapper.h>

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

struct WorkspaceMethodInternalRecord {
  std::string desc;
  std::vector<std::string> author;
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

  [[nodiscard]] int count_overloads() const;
  [[nodiscard]] std::vector<std::vector<std::string>> generic_overloads() const;
  [[nodiscard]] bool has_any() const;
  [[nodiscard]] bool has_overloads() const;
  [[nodiscard]] std::string docstring() const;
  [[nodiscard]] std::string header(const std::string& name, int = 0) const;
  [[nodiscard]] std::string call(const std::string& name) const;
};

const std::unordered_map<std::string, WorkspaceMethodInternalRecord>&
internal_workspace_methods();
