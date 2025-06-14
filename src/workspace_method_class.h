#pragma once

#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "format_tags.h"
#include <auto_wsv_value_wrapper.h>

inline constexpr char named_input_prefix = '@';
inline constexpr char internal_prefix    = '_';

class Method {
  std::string name{};
  std::vector<std::string> outargs{};
  std::vector<std::string> inargs{};
  std::optional<Wsv> setval{std::nullopt};
  bool overwrite_setval{false};

 public:
  Method();
  Method(const std::string& name,
         const std::vector<std::string>& args,
         const std::unordered_map<std::string, std::string>& kwargs);
  Method(std::string name, const Wsv& wsv, bool = false);
  Method(const std::string& name,
         const std::vector<std::string>& ins,
         const std::vector<std::string>& outs,
         const std::optional<Wsv>& wsv,
         bool overwrite);

  [[nodiscard]] const std::string& get_name() const { return name; }
  [[nodiscard]] const std::vector<std::string>& get_outs() const {
    return outargs;
  }
  [[nodiscard]] const std::vector<std::string>& get_ins() const {
    return inargs;
  }
  [[nodiscard]] const std::optional<Wsv>& get_setval() const { return setval; }
  [[nodiscard]] bool overwrite() const { return overwrite_setval; }

  void operator()(Workspace& ws) const;
  void add_defaults_to_agenda(Agenda& agenda) const;

  [[nodiscard]] std::string sphinx_list_item() const;

  friend std::ostream& operator<<(std::ostream& os, const Method& m);
};

template <>
struct std::formatter<Wsv> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  [[nodiscard]] std::string to_string(const Wsv& wsv) const;

  template <class FmtContext>
  FmtContext::iterator format(const Wsv& v, FmtContext& ctx) const {
    if (tags.short_str) {
      return std::format_to(
          ctx.out(), "{}{}{}", tags.quote(), v.type_name(), tags.quote());
    }

    return tags.format(ctx, to_string(v));
  }
};

template <>
struct std::formatter<Method> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Method& v, FmtContext& ctx) const {
    tags.add_if_bracket(ctx, '[');
    tags.add_if_bracket(ctx, '\n');

    const std::string_view sep   = tags.sep(true);
    const std::string_view quote = tags.quote();

    if (tags.names) {
      const std::string_view ow = v.overwrite() ? " (overwrite)"sv : ""sv;
      tags.format(ctx, quote, v.get_name(), ow, quote, sep);
    }

    if (v.get_setval().has_value()) {
      tags.format(ctx, v.get_setval().value());
    } else {
      tags.format(ctx,
                  quote,
                  "outputs"sv,
                  quote,
                  ": "sv,
                  v.get_outs(),
                  sep,
                  quote,
                  "inputs"sv,
                  quote,
                  ": "sv,
                  v.get_ins());
    }

    tags.add_if_bracket(ctx, '\n');
    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};

template <>
struct std::formatter<Agenda> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Agenda& v, FmtContext& ctx) const {
    tags.add_if_bracket(ctx, '[');

    const std::string_view sep   = tags.sep(true);
    const std::string_view quote = tags.quote();

    if (tags.names) {
      tags.format(ctx,
                  quote,
                  v.get_name(),
                  quote,
                  " (checked: "sv,
                  v.is_checked(),
                  ')');
    }

    tags.format(ctx,
                sep,
                v.get_methods(),
                sep,
                "  Shared: "sv,
                v.get_share(),
                sep,
                "  Copied: "sv,
                v.get_copy());

    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};
