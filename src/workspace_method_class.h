#pragma once

#include <optional>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "auto_wsg.h"
#include "format_tags.h"

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

  friend std::ostream& operator<<(std::ostream& os, const Method& m);
};

template <>
struct std::formatter<Wsv> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Wsv& v, FmtContext& ctx) const {
    if (tags.short_str) {
      return std::format_to(
          ctx.out(), "{}{}{}", tags.quote(), v.type_name(), tags.quote());
    }

    std::visit(
        [*this, &ctx]<typename T>(const std::shared_ptr<T>& x) {
          std::formatter<T> fmt{};
          this->make_compat(fmt);
          fmt.format(*x, ctx);
        },
        v.value);

    return ctx.out();
  }
};

template <>
struct std::formatter<Method> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Method& v, FmtContext& ctx) const {
    if (tags.bracket) std::format_to(ctx.out(), "[\n"sv);

    const std::string_view sep = tags.sep(true);

    if (tags.names) {
      const std::string_view ow = v.overwrite() ? " (overwrite)"sv : ""sv;
      std::format_to(ctx.out(),
                     "{}{}{}{}{}",
                     tags.quote(),
                     v.get_name(),
                     tags.quote(),
                     ow,
                     sep);
    }

    if (v.get_setval().has_value()) {
      std::formatter<Wsv> setval{};
      make_compat(setval);
      setval.format(v.get_setval().value(), ctx);
    } else {
      std::formatter<std::vector<std::string>> outargs{};
      std::formatter<std::vector<std::string>> inargs{};
      make_compat(outargs, inargs);
      std::format_to(ctx.out(), "outputs: "sv);
      outargs.format(v.get_outs(), ctx);
      std::format_to(ctx.out(), "{}inputs: ", sep);
      inargs.format(v.get_ins(), ctx);
    }

    if (tags.bracket) std::format_to(ctx.out(), "\n]"sv);
    return ctx.out();
  }
};

template <>
struct std::formatter<Agenda> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Agenda& v, FmtContext& ctx) const {
    tags.add_if_bracket(ctx, '[');

    const std::string_view sep = tags.sep(true);

    if (tags.names) {
      std::format_to(ctx.out(),
                     "{}{}{} (checked: {}):",
                     tags.quote(),
                     v.get_name(),
                     tags.quote(),
                     v.is_checked());
    }

    std::formatter<std::vector<Method>> methods{};
    std::formatter<std::vector<std::string>> share{};
    std::formatter<std::vector<std::string>> copy{};
    make_compat(methods, share, copy);

    methods.format(v.get_methods(), ctx);
    std::format_to(ctx.out(), "{}  Shared: ", sep);
    share.format(v.get_share(), ctx);
    std::format_to(ctx.out(), "{}  Copied: ", sep);
    copy.format(v.get_copy(), ctx);

    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};
