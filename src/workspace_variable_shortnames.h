#pragma once

#include <string>
#include <unordered_map>

struct WsvShortForm {
  std::string desc;
  std::string title{"Normal Variables"};
};

const std::unordered_map<std::string, WsvShortForm>&
workspace_variables_shortnames();
