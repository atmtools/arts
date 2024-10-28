#pragma once

#include <string>
#include <unordered_map>

struct WorkspaceGroupRecord {
  std::string file;
  std::string desc;

  //! Set to level of array nesting (0 is not an array)
  int array_depth{0};

  //! Set to true if python must treat this as a value
  bool value_type{false};

  //! Set to true if the type is a map of some type
  bool map_type{false};
};

const std::unordered_map<std::string, WorkspaceGroupRecord>&
internal_workspace_groups();
