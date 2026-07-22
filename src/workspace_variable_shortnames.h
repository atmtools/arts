#pragma once

#include <format_tags.h>

#include <string>
#include <unordered_map>

struct WsvShortForm {
  std::string desc;
  std::string title{"Normal Variables"};
};

const std::unordered_map<std::string, WsvShortForm>& workspace_variables_shortnames();

template <> struct std::formatter<WsvShortForm> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const WsvShortForm& wsg, FmtContext& ctx) const {
    return tags.format(ctx, "WsvShortForm{\n  .desc="sv, wsg.desc, "\n  .title="sv, wsg.title, "\n}"sv);
  }
};
