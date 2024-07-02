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

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Wsv& v, FmtContext& ctx) const {
    if (tags.short_str) {
      std::format_to(ctx, "{}", v.type_name());
    } else {
      std::visit(
          [*this, &ctx]<typename T>(const std::shared_ptr<T>& x) {
            if constexpr (arts_formattable<T>) {
              std::formatter<T> fmt{};
              make_compat(fmt);
              fmt.format(*x, ctx);
            } else {
              std::ranges::copy("Not formattable"sv, ctx.out());
            }
          },
          v.value);
    }
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

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Method& v, FmtContext& ctx) const {
    if (tags.bracket) std::ranges::copy("[\n"sv, ctx.out());

    const std::string_view sep = tags.comma ? ",\n"sv : "\n"sv;

    if (tags.names) {
      const std::string_view quote = tags.bracket ? R"(")" : ""sv;
      std::format_to(ctx, "{}{}{}", quote, v.get_name(), quote);
      if (v.overwrite()) std::format_to(ctx, " (overwrite)");
      std::format_to(ctx, "{}", sep);
    }

    if (v.get_setval().has_value()) {
      std::formatter<Wsv> setval{};
      make_compat(setval);
      setval.format(v.get_setval().value(), ctx);
    } else {
      std::formatter<std::vector<std::string>> outargs{};
      std::formatter<std::vector<std::string>> inargs{};
      make_compat(outargs, inargs);
      std::format_to(ctx, "{}outputs: ", sep);
      outargs.format(v.get_outs(), ctx);
      std::format_to(ctx, "{}inputs: ", sep);
      inargs.format(v.get_ins(), ctx);
    }

    if (tags.bracket) std::ranges::copy("\n]"sv, ctx.out());
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

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Agenda& v, FmtContext& ctx) const {
    if (tags.bracket) std::ranges::copy("["sv, ctx.out());

    const std::string_view sep = tags.comma ? ",\n"sv : "\n"sv;

    if (tags.names) {
      const std::string_view quote = tags.bracket ? R"(")" : ""sv;
      std::format_to(ctx,
                     "{}{}{} (checked: {}):\n",
                     quote,
                     v.get_name(),
                     quote,
                     v.is_checked(),
                     sep);
    }

    std::formatter<std::vector<Method>> methods{};
    std::formatter<std::vector<std::string>> share{};
    std::formatter<std::vector<std::string>> copy{};
    make_compat(methods, share, copy);

    methods.format(v.get_methods(), ctx);
    std::format_to(ctx, "{}Shared: ", sep);
    share.format(v.get_methods(), ctx);
    std::format_to(ctx, "{}Copied: ", sep);
    copy.format(v.get_methods(), ctx);

    if (tags.bracket) std::ranges::copy("]"sv, ctx.out());
    return ctx.out();
  }
};
