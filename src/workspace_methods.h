#pragma once

#include <wsv_value_wrapper.h>

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

struct WorkspaceMethodInternalRecord {
  std::string                     desc;
  std::vector<std::string>        author;
  std::string                     return_type{"void"};
  std::string                     return_desc{};
  std::vector<std::string>        out{};
  std::vector<std::string>        gout{};
  std::vector<std::string>        gout_type{};
  std::vector<std::string>        gout_desc{};
  std::vector<std::string>        in{};
  std::vector<std::string>        gin{};
  std::vector<std::string>        gin_type{};
  std::vector<std::optional<Wsv>> gin_value{};
  std::vector<std::string>        gin_desc{};
  bool                            pass_workspace{false};

  [[nodiscard]] int                                   count_overloads() const;
  [[nodiscard]] std::vector<std::vector<std::string>> generic_overloads() const;
  [[nodiscard]] bool                                  has_any() const;
  [[nodiscard]] bool                                  has_overloads() const;
  [[nodiscard]] std::string                           docstring() const;
  [[nodiscard]] std::string                           header(const std::string& name, int = 0) const;
  [[nodiscard]] std::string                           call(const std::string& name) const;
};

const std::unordered_map<std::string, WorkspaceMethodInternalRecord>& internal_workspace_methods();

template <> struct std::formatter<WorkspaceMethodInternalRecord> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const WorkspaceMethodInternalRecord& wsm, FmtContext& ctx) const {
    return tags.format(ctx,
                       "WorkspaceMethodInternalRecord{"sv,
                       "\n  .desc="sv,
                       wsm.desc,
                       "\n  .author="sv,
                       wsm.author,
                       "\n  .return_type="sv,
                       wsm.return_type,
                       "\n  .return_desc="sv,
                       wsm.return_desc,
                       "\n  .out="sv,
                       wsm.out,
                       "\n  .gout="sv,
                       wsm.gout,
                       "\n  .gout_type="sv,
                       wsm.gout_type,
                       "\n  .gout_desc="sv,
                       wsm.gout_desc,
                       "\n  .in="sv,
                       wsm.in,
                       "\n  .gin="sv,
                       wsm.gin,
                       "\n  .gin_type="sv,
                       wsm.gin_type,
                       "\n  .gin_value="sv,
                       wsm.gin_value,
                       "\n  .gin_desc="sv,
                       wsm.gin_desc,
                       "\n  .pass_workspace="sv,
                       wsm.pass_workspace ? "true"sv : "false"sv,
                       "\n}"sv);
  }
};
