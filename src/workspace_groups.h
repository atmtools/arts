#pragma once

#include <format_tags.h>

#include <string>
#include <unordered_map>
#include <vector>

struct WorkspaceGroupRecord {
  std::string file;
  std::string desc;

  //! Set to true if python must treat this as a value
  bool value_type{false};

  //! Set to true if the type is a map of some type
  bool map_type{false};
};

const std::unordered_map<std::string, WorkspaceGroupRecord>& internal_workspace_groups();

void add_arrays_of(std::unordered_map<std::string, WorkspaceGroupRecord>& wsg_data,
                   const std::vector<std::string>&                        types,
                   std::vector<std::string>                               extra_headers);

template <> struct std::formatter<WorkspaceGroupRecord> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const WorkspaceGroupRecord& wsg, FmtContext& ctx) const {
    return tags.format(ctx,
                       "WorkspaceGroupRecord{"sv,
                       "\n  .file="sv,
                       wsg.file,
                       "\n  .desc="sv,
                       wsg.desc,
                       "\n  .value_type="sv,
                       wsg.value_type ? "true"sv : "false"sv,
                       "\n  .map_type="sv,
                       wsg.map_type ? "true"sv : "false"sv,
                       "\n}"sv);
  }
};
