#ifndef callback_h
#define callback_h

class Workspace;

#include <format_tags.h>

#include <functional>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

struct CallbackOperator {
  std::function<void(const std::shared_ptr<Workspace>&)> callback{};
  std::vector<std::string> inputs{};
  std::vector<std::string> outputs{};

  void operator()(Workspace& ws) const;

  friend std::ostream& operator<<(std::ostream& os, const CallbackOperator&);
};

template <>
struct std::formatter<CallbackOperator> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x._make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const CallbackOperator& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.comma ? ", "sv : " "sv;

    std::formatter<std::vector<std::string>> strs;

    make_compat(strs);

    if (tags.bracket) std::ranges::copy("["sv, ctx.out());
    std::format_to(ctx, "input: ");
    strs.format(v.inputs, ctx);
    std::format_to(ctx, "{}output: ", sep);
    strs.format(v.outputs, ctx);
    if (tags.bracket) std::ranges::copy("]"sv, ctx.out());

    return ctx.out();
  }
};

#endif
