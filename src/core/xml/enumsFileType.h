#pragma once

#include "enums-common-helper.h"

enum class FileType : unsigned char {
  ascii,
  zascii,
  binary,
};

template <> constexpr bool good_enum<FileType>(FileType x) noexcept {
  const auto v = static_cast<std::size_t>(x);
  return v < 3;
}

template <> struct enumdocs<FileType> {
  static std::string_view           str() noexcept;
  static constexpr std::string_view name = "FileType"sv;
};

namespace enumtyps {
inline constexpr std::array FileTypeTypes = {
    FileType::ascii,
    FileType::zascii,
    FileType::binary,
};
}  // namespace enumtyps

namespace enumstrs {

template <> struct enum_str_data<FileType, 0> {
  static constexpr std::array strs = {
      "ascii"sv,
      "zascii"sv,
      "binary"sv,
  };
};

template <> struct enum_str_data<FileType, 1> {
  static constexpr std::array strs = {
      "ASCII"sv,
      "ZASCII"sv,
      "BINARY"sv,
  };
};

template <> struct enum_str_data<FileType, 2> {
  static constexpr std::array strs = {
      "Ascii"sv,
      "Zip"sv,
      "Binary"sv,
  };
};

template <> struct enum_str_data<FileType, 3> {
  static constexpr std::array strs = {
      "text"sv,
      "zip"sv,
      "bin"sv,
  };
};

template <int i = 0> inline constexpr auto FileTypeNames = enum_str_data<FileType, i>::strs;
}  // namespace enumstrs

template <int i = 0> constexpr std::string_view toString(FileType x) requires(i >= 0 and i < 4) {
  if (good_enum(x)) return enumstrs::FileTypeNames<i>[static_cast<std::size_t>(x)];
  return "BAD FileType"sv;
}

template <> constexpr FileType to<FileType>(const std::string_view x) {
  switch (enumhash(x)) {
    case enumhash("ascii"sv):
      if (x == "ascii"sv) return FileType::ascii;
      break;
    case enumhash("zascii"sv):
      if (x == "zascii"sv) return FileType::zascii;
      break;
    case enumhash("binary"sv):
      if (x == "binary"sv) return FileType::binary;
      break;
    case enumhash("ASCII"sv):
      if (x == "ASCII"sv) return FileType::ascii;
      break;
    case enumhash("ZASCII"sv):
      if (x == "ZASCII"sv) return FileType::zascii;
      break;
    case enumhash("BINARY"sv):
      if (x == "BINARY"sv) return FileType::binary;
      break;
    case enumhash("Ascii"sv):
      if (x == "Ascii"sv) return FileType::ascii;
      break;
    case enumhash("Zip"sv):
      if (x == "Zip"sv) return FileType::zascii;
      break;
    case enumhash("Binary"sv):
      if (x == "Binary"sv) return FileType::binary;
      break;
    case enumhash("text"sv):
      if (x == "text"sv) return FileType::ascii;
      break;
    case enumhash("zip"sv):
      if (x == "zip"sv) return FileType::zascii;
      break;
    case enumhash("bin"sv):
      if (x == "bin"sv) return FileType::binary;
      break;
  }
  throw std::runtime_error(arts::format(R"-x-(Bad input "{}"

See https://atmtools.github.io/arts-docs-master/pyarts3.arts.FileType.html for valid options.
)-x-",
                                        x));
}

namespace enumsize {
inline constexpr std::size_t FileTypeSize = 3;
}

template <> struct std::formatter<FileType> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const FileType& v, FmtContext& ctx) const {
    return tags.format(ctx, tags.quote(), toString<0>(v), tags.quote());
  }
};
