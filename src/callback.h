#ifndef callback_h
#define callback_h

class Workspace;

#include <format_tags.h>

#include <functional>
#include <string>
#include <string_view>
#include <vector>

struct CallbackOperator {
  std::function<void(Workspace&)> callback{};
  std::vector<std::string> inputs{};
  std::vector<std::string> outputs{};

  void operator()(Workspace& ws) const;
};

template <>
struct std::formatter<CallbackOperator> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const CallbackOperator& v,
                              FmtContext& ctx) const {
    const std::string_view quote = tags.quote();

    tags.add_if_bracket(ctx, '{');

    tags.format(ctx,
                quote,
                "input"sv,
                quote,
                ": "sv,
                v.inputs,
                tags.sep(true),
                quote,
                "output"sv,
                quote,
                ": "sv,
                v.outputs);

    tags.add_if_bracket(ctx, '}');
    return ctx.out();
  }
};

template <>
struct xml_io_stream<CallbackOperator> {
  static constexpr std::string_view type_name = "CallbackOperator"sv;

  static void write(std::ostream& os,
                    const CallbackOperator& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   CallbackOperator& x,
                   bifstream* pbifs = nullptr);
};

#endif
