#pragma once

#include <string>
#include <unordered_map>

struct WsvShortForm {
  std::string related_to;
  std::string desc;

  // Variable naming style.  Used for table titles in documentation.
  std::string style{"Normal Variable"};
};

const std::unordered_map<std::string, WsvShortForm>&
workspace_variables_shortnames();
