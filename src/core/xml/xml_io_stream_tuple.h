#pragma once

#include <xml_io_base.h>
#include <xml_io_stream.h>

#include <concepts>
#include <tuple>

template <typename... Ts>
struct xml_io_stream_name<std::tuple<Ts...>> {
  static constexpr std::string_view name = "Tuple"sv;
};

template <arts_xml_ioable... Ts>
struct xml_io_stream<std::tuple<Ts...>> {
  constexpr static std::string_view type_name =
      xml_io_stream_name_v<std::tuple<Ts...>>;

  static void write(std::ostream& os,
                    const std::tuple<Ts...>& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) {
    std::println(
        os, R"(<{0} name="{1}" ntypes="{2}">)", type_name, name, sizeof...(Ts));

    std::apply(
        [&](const Ts&... t) { (xml_io_stream<Ts>::write(os, t, pbofs), ...); },
        n);

    std::println(os, R"(</{0}>)", type_name);
  }

  static void read(std::istream& is,
                   std::tuple<Ts...>& n,
                   bifstream* pbifs = nullptr) {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);

    Size ntypes;
    tag.get_attribute_value("ntypes", ntypes);

    if (ntypes != sizeof...(Ts))
      throw std::runtime_error(std::format(
          "Expecting only {} items, got {}", sizeof...(Ts), ntypes));

    std::apply([&](Ts&... t) { (xml_io_stream<Ts>::read(is, t, pbifs), ...); },
               n);

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  }
};

template <typename A, typename B>
struct xml_io_stream_name<std::pair<A, B>> {
  static constexpr std::string_view name = "Pair"sv;
};

template <arts_xml_ioable A, arts_xml_ioable B>
struct xml_io_stream<std::pair<A, B>> {
  constexpr static std::string_view type_name =
      xml_io_stream_name_v<std::pair<A, B>>;

  static void write(std::ostream& os,
                    const std::pair<A, B>& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) {
    std::println(os, R"(<{0} name="{1}" ntypes="{2}">)", type_name, name, 2);

    xml_io_stream<A>::write(os, n.first, pbofs, "first");
    xml_io_stream<A>::write(os, n.second, pbofs, "second");

    std::println(os, R"(</{0}>)", type_name);
  }

  static void read(std::istream& is,
                   std::pair<A, B>& n,
                   bifstream* pbifs = nullptr) {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);

    Size ntypes;
    tag.get_attribute_value("ntypes", ntypes);

    if (ntypes != 2)
      throw std::runtime_error(
          std::format("Expecting only {} items, got {}", 2, ntypes));

    xml_io_stream<A>::read(is, n.first, pbifs);
    xml_io_stream<A>::read(is, n.second, pbifs);

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  }
};
