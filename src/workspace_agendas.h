#pragma once

#include <string>
#include <unordered_map>
#include <vector>

struct StringVectorAgendaHelper {
  std::string test{};
  std::string constraint{};
  std::vector<std::string> printables{};

  constexpr StringVectorAgendaHelper() = default;

  template <typename... Args>
  constexpr StringVectorAgendaHelper(std::string&& Test,
                                     std::string&& Constraint,
                                     Args&&... Printables)
      : test(std::move(Test)),
        constraint(std::move(Constraint)),
        printables(
            std::vector<std::string>{std::forward<Args>(Printables)...}) {}
};

struct WorkspaceAgendaInternalRecord {
  std::string desc{};
  std::vector<std::string> output{};
  std::vector<std::string> input{};
  std::vector<std::string> enum_options{};
  std::string enum_default{};
  std::vector<StringVectorAgendaHelper> output_constraints{};
};

const std::unordered_map<std::string, WorkspaceAgendaInternalRecord>&
internal_workspace_agendas();
