#include "workspace_agendas.h"

std::unordered_map<std::string, WorkspaceAgendaInternalRecord> internal_workspace_agendas() {
  std::unordered_map<std::string, WorkspaceAgendaInternalRecord> wsa_data;

  wsa_data["y_out_agenda"] = {
    .desc = R"--(PURE OUT
)--",
    .output = {"y"}
    };
    
  wsa_data["x_out_agenda"] = {
    .desc = R"--(IN / OUT
)--",
    .output = {"x"},
    .input = {"y"}
    };

  return wsa_data;
}
