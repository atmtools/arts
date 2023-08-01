#pragma once

#include <string>
#include <unordered_map>
#include <vector>

struct WorkspaceAgendaInternalRecord {
  std::string desc;
  std::vector <std::string> output;
  std::vector <std::string> input;
  bool array{false};
};

std::unordered_map<std::string, WorkspaceAgendaInternalRecord> internal_workspace_agendas();
