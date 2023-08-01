#include "workspace_variables.h"

#include <iomanip>

std::unordered_map<std::string, WorkspaceVariableInternalRecord> internal_workspace_variables() {
  std::unordered_map<std::string, WorkspaceVariableInternalRecord> wsv_data;
  
  wsv_data["y"] = {
    .desc          = R"--(A variable
)--",
    .type          = "Vector",
    .default_value = std::nullopt
  };
  
  wsv_data["x"] = {
    .desc          = R"--(Another variable
)--",
    .type          = "Numeric",
    .default_value = Numeric{1.0}
  };
  
  wsv_data["num"] = {
    .desc          = R"--(A third variable
)--",
    .type          = "Index"
  };

  return wsv_data;
}
