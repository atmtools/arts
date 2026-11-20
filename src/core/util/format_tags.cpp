#include "format_tags.h"

#include <cstdio>
#include <ostream>

std::string format_tags::get_format_args() const {
  std::string buf{"{:"};
  buf.reserve(16);

  if (names) buf.push_back('N');
  if (bracket) buf.push_back('B');
  if (quoted) buf.push_back('q');
  if (short_str) buf.push_back('s');
  if (comma) buf.push_back(',');
  if (newline) buf.push_back('n');
  if (help) buf.push_back('h');
  if (io) buf += "IO"sv;

  buf.push_back('}');
  if (buf.empty()) return buf;

  return buf;
}

std::string_view format_tags::sep() const {
  if (newline) {
    if (comma) return ",\n"sv;
    return "\n"sv;
  }
  if (comma) return ", "sv;
  return " "sv;
}

std::string_view format_tags::quote() const {
  if (quoted) return R"(")"sv;
  return ""sv;
}

template <> void format_tags::add_if_bracket(std::format_context& ctx, std::string_view x) const try {
  if (bracket) std::format_to(ctx.out(), "{}", x);
} catch (const std::exception& e) {
  throw std::runtime_error("Error in single_format with fmt-string: " + get_format_args() + "\n" + e.what());
}

template <> std::string format_tags::vformat(const std::string& x) const { return x; }

template <> std::string format_tags::vformat(const std::string_view& x) const { return std::string(x); }

namespace arts {
std::string vformat(std::string_view fmt, std::format_args args) { return std::vformat(fmt, args); }

std::format_context::iterator vformat_to(std::format_context::iterator out,
                                         std::string_view              fmt,
                                         std::format_args              args) {
  return std::vformat_to(out, fmt, args);
}

void vprint(std::FILE* stream, std::string_view fmt, std::format_args args, bool newline) {
  if (newline)
    std::println(stream, "{}", std::vformat(fmt, args));
  else
    std::vprint_unicode(stream, fmt, args);
}

void vprint(std::ostream& os, std::string_view fmt, std::format_args args, bool newline) {
  if (newline)
    std::println(os, "{}", std::vformat(fmt, args));
  else
    std::vprint_unicode(os, fmt, args);
}

void println(std::FILE* stream) { std::print(stream, "\n"); }

void println(std::ostream& os) { std::print(os, "\n"); }
}  // namespace arts

template <> void format_tags::single_format(std::format_context& ctx, const Numeric& x) const {
  std::formatter<Numeric>{}.format(x, ctx);
}

template <> void format_tags::single_format(std::format_context& ctx, const Index& x) const {
  std::formatter<Index>{}.format(x, ctx);
}

template <> void format_tags::single_format(std::format_context& ctx, const Size& x) const {
  std::formatter<Size>{}.format(x, ctx);
}

template <> void format_tags::single_format(std::format_context& ctx, const int& x) const {
  std::formatter<int>{}.format(x, ctx);
}

template <> void format_tags::single_format(std::format_context& ctx, const bool& x) const {
  std::formatter<bool>{}.format(x, ctx);
}

template <> void format_tags::single_format(std::format_context& ctx, const std::string& x) const {
  std::formatter<std::string>{}.format(x, ctx);
}

template <> void format_tags::single_format(std::format_context& ctx, const std::string_view& x) const {
  std::formatter<std::string_view>{}.format(x, ctx);
}
