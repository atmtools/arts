#pragma once

#include <string>
#include <unordered_map>

struct WorkspaceGroupRecord {
  std::string file;
  std::string desc;

  //! Set to true if python must treat this as a value
  bool value_type{false};

  //! Set to true if this is an array of values
  bool is_array{false};

  //! Set to true if this is a simple contiguous array (simple means of one of the basic types of C++, i.e., double or int or similiar)
  bool is_simple_contiguous{false};
};

const std::unordered_map<std::string, WorkspaceGroupRecord>&
internal_workspace_groups();
