#pragma once

#include <aggregate_helper.h>

#include "xml_io_base.h"
#include "xml_io_stream.h"

template <typename T> struct xml_io_stream_aggregate {
  constexpr static bool value = false;
};

template <typename T> constexpr bool xml_io_stream_aggregate_v = xml_io_stream_aggregate<T>::value;

template <typename T>
concept xml_io_aggregratable = xml_io_stream_aggregate_v<T>;

template <xml_io_aggregratable T> struct xml_io_stream<T> {
  constexpr static std::string_view type_name = xml_io_stream_name_v<T>;

  static_assert(arts_aggregate<T>, "xml_io_aggregratable types must be arts_aggregate");

  using _lambda = decltype([](auto& v) { return as_tuple(v); });
  using ctup_t  = std::invoke_result_t<_lambda, const T&>;
  using mtup_t  = std::invoke_result_t<_lambda, T&>;
  using inner   = xml_io_stream<mtup_t>;

  static void write(std::ostream& os, const T& t, bofstream* pbofs = nullptr, std::string_view name = ""sv) {
    XMLTag tag(type_name, "name", name);
    tag.write_to_stream(os);

    std::apply([&os, &pbofs]<typename... Ts>(const Ts&... v) { (xml_io_stream<Ts>::write(os, v, pbofs), ...); },
               as_tuple(t));

    tag.write_to_end_stream(os);
  }

  static void read(std::istream& is, T& t, bifstream* pbifs = nullptr) try {
    XMLTag tag{};
    tag.read_from_stream(is);
    tag.check_name(type_name);

    std::apply([&is, &pbifs]<typename... Ts>(Ts&... v) { (xml_io_stream<Ts>::read(is, v, pbifs), ...); }, as_tuple(t));

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::format("Cannot read {}:\n{}", type_name, e.what()));
  }
};
