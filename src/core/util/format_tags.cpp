#include "format_tags.h"

std::string format_tags::get_format_args() const {
  std::string buf{"{:"};
  buf.reserve(16);

  if (names)     buf.push_back('N');
  if (bracket)   buf.push_back('B');
  if (quoted)    buf.push_back('q');
  if (short_str) buf.push_back('s');
  if (comma)     buf.push_back(',');
  
  if (io)        buf += "IO"sv;

  buf.push_back('}');
  if (buf.empty()) return buf;

  return buf;
}

std::string_view format_tags::sep(bool newline) const {
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
