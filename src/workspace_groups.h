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

  //! Set to true if this is a simple contiguous array (simple means of one of the basic types of C++, i.e., double or int or similiar)
  bool is_simple_contiguous{false};

  //! Set to true if this cannot be automatically initialized in python
  bool skip_pyinit{false};
};

const std::unordered_map<std::string, WorkspaceGroupRecord>&
internal_workspace_groups();
