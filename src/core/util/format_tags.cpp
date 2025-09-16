#include "format_tags.h"

std::string format_tags::get_format_args() const {
  std::string buf{"{:"};
  buf.reserve(16);

  if (names) buf.push_back('N');
  if (bracket) buf.push_back('B');
  if (quoted) buf.push_back('q');
  if (short_str) buf.push_back('s');
  if (comma) buf.push_back(',');
  if (newline) buf.push_back('n');
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

template <>
void format_tags::add_if_bracket(std::format_context& ctx, char x) const try {
  if (bracket) std::format_to(ctx.out(), "{}", x);
} catch (const std::exception& e) {
  throw std::runtime_error("Error in single_format with fmt-string: " +
                           get_format_args() + "\n" + e.what());
}

template <>
std::string format_tags::vformat(const std::string& x) const {
  return x;
}

template <>
std::string format_tags::vformat(const std::string_view& x) const {
  return std::string(x);
}
