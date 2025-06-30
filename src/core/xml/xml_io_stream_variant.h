#pragma once

#include <xml_io_base.h>
#include <xml_io_stream.h>

#include <concepts>
#include <variant>

//! Try overloading these if it ever gets too slow.
namespace {
template <typename... Ts>
constexpr static std::array type_names{xml_io_stream<Ts>::type_name...};

template <Size N>
consteval bool unique_names(std::array<std::string_view, N> x) {
  std::ranges::sort(x);
  return std::adjacent_find(x.begin(), x.end()) == x.end();
}

template <typename... Ts>
constexpr std::string_view variant_unique_name(const std::variant<Ts...>& x) {
  return type_names<Ts...>[x.index()];
}

template <Size I = 0, typename... Ts>
std::variant<Ts...> variant_type_init(const std::string_view& x) {
  if (type_names<Ts...>[I] == x) {
    return std::variant_alternative_t<I, std::variant<Ts...>>{};
  }

  if constexpr (I + 1 < sizeof...(Ts)) {
    return variant_type_init<I + 1, Ts...>(x);
  } else {
    throw std::runtime_error(
        std::format(R"(Cannot understand the variant type: "{}")", x));
  }
}

template <typename... Ts>
bool variant_write(std::ostream& os,
                   const std::variant<Ts...>& n,
                   bofstream* pbofs) {
  const auto call = []<typename T>(std::ostream& os,
                                   const T* const e,
                                   bofstream* pbofs) -> bool {
    if (e) xml_io_stream<T>::write(os, *e, pbofs);
    return e;
  };
  return (call(os, std::get_if<Ts>(&n), pbofs) or ...);
}

template <typename... Ts>
bool variant_read(std::istream& is, std::variant<Ts...>& n, bifstream* pbifs) {
  const auto call = []<typename T>(
                        std::istream& is, T* e, bifstream* pbifs) -> bool {
    if (e) xml_io_stream<T>::read(is, *e, pbifs);
    return e;
  };
  return (call(is, std::get_if<Ts>(&n), pbifs) or ...);
}

template <typename... Ts>
concept uniquely_variant = unique_names(type_names<Ts...>);
}  // namespace

template <arts_xml_ioable... Ts>
  requires(uniquely_variant<Ts...>)
struct xml_io_stream<std::variant<Ts...>> {
  constexpr static std::string_view type_name = "Variant"sv;

  static void write(std::ostream& os,
                    const std::variant<Ts...>& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) {
    std::println(os,
                 R"(<{0} type="{1}" name="{2}">)",
                 type_name,
                 variant_unique_name(n),
                 name);

    if (not variant_write(os, n, pbofs)) {
      throw std::runtime_error("Failed to write, got nullptr");
    }

    std::println(os, R"(</{0}>)", type_name);
  }

  static void read(std::istream& is,
                   std::variant<Ts...>& n,
                   bifstream* pbifs = nullptr) {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);

    String type;
    tag.get_attribute_value("type", type);

    if (not variant_read(is, n = variant_type_init<0, Ts...>(type), pbifs)) {
      throw std::runtime_error("Failed to read, got nullptr");
    }

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  }
};
