#pragma once

#include <format_tags.h>

#include <string>
#include <unordered_map>
#include <vector>

struct StringVectorAgendaHelper {
  std::string              test{};
  std::string              constraint{};
  std::vector<std::string> printables{};

  constexpr StringVectorAgendaHelper() = default;

  template <typename... Args>
  constexpr StringVectorAgendaHelper(std::string&& Test, std::string&& Constraint, Args&&... Printables)
      : test(std::move(Test)),
        constraint(std::move(Constraint)),
        printables(std::vector<std::string>{std::forward<Args>(Printables)...}) {}
};

struct WorkspaceAgendaInternalRecord {
  std::string                           desc{};
  std::vector<std::string>              output{};
  std::vector<std::string>              input{};
  std::vector<std::string>              enum_options{};
  std::string                           enum_default{};
  std::vector<StringVectorAgendaHelper> output_constraints{};
  std::string                           named_operator{};
};

const std::unordered_map<std::string, WorkspaceAgendaInternalRecord>& internal_workspace_agendas();

template <> struct std::formatter<StringVectorAgendaHelper> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const StringVectorAgendaHelper& wair, FmtContext& ctx) const {
    return tags.format(ctx,
                       "StringVectorAgendaHelper{\n  .test="sv,
                       wair.test,
                       "\n  .constraint="sv,
                       wair.constraint,
                       "\n  .printables="sv,
                       wair.printables,
                       "\n}"sv);
  }
};

template <> struct std::formatter<WorkspaceAgendaInternalRecord> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const WorkspaceAgendaInternalRecord& wair, FmtContext& ctx) const {
    return tags.format(ctx,
                       "WorkspaceAgendaInternalRecord{\n  .desc="sv,
                       wair.desc,
                       "\n  .output="sv,
                       wair.output,
                       "\n  .input="sv,
                       wair.input,
                       "\n  .enum_options="sv,
                       wair.enum_options,
                       "\n  .enum_default="sv,
                       wair.enum_default,
                       "\n  .output_constraints="sv,
                       wair.output_constraints,
                       "\n  .named_operator="sv,
                       wair.named_operator,
                       "\n}"sv);
  }
};
