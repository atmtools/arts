#include "workspace_methods.h"

#include <limits>
#include <optional>

std::unordered_map<std::string, WorkspaceMethodInternalRecord>
internal_workspace_methods() {
  std::unordered_map<std::string, WorkspaceMethodInternalRecord> wsm_data;

  wsm_data["print"] = {
    .desc = R"--(PURE GIN
)--",
    .author = {"Richard Larsson"},
    .gin = {"v"},
    .gin_type = {"Any"},
    .gin_value = {std::nullopt},
    .gin_desc = {"A variable"},
  };

  wsm_data["create"] = {
    .desc = R"--(PURE OUT
)--",
    .author = {"Richard Larsson"},
    .out = {"num"}
  };

  wsm_data["add"] = {
    .desc = R"--(INOUT
)--",
    .author = {"Richard Larsson"},
    .out = {"y"},
    .in = {"y", "x"}
  };

  wsm_data["new_y"] = {
    .desc = R"--(OUT / IN
)--",
    .author = {"Richard Larsson"},
    .out = {"y"},
    .in = {"x"}
  };

  wsm_data["append"] = {
    .desc = R"--(GOUT / GIN --- GENERIC
)--",
    .author = {"Richard Larsson"},
    .gout = {"v_out"},
    .gout_type = {"Vector"},
    .gout_desc = {"A variable"},
    .gin = {"v_in", "n"},
    .gin_type = {"Vector", "Index, Numeric"},
    .gin_value = {std::nullopt, std::nullopt},
    .gin_desc = {"Vector input", "Something to append"},
  };

  wsm_data["defval_append"] = {
    .desc = R"--(GOUT / GIN --- GENERIC
)--",
    .author = {"Richard Larsson"},
    .gout = {"v_out"},
    .gout_type = {"Vector"},
    .gout_desc = {"A variable"},
    .gin = {"v_in", "n"},
    .gin_type = {"Vector", "Index, Numeric"},
    .gin_value = {std::nullopt, Index{1}},
    .gin_desc = {"Vector input", "Something to append"},
  };

  return wsm_data;
}