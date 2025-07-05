#pragma once

#include <concepts>
#include <tuple>

#include "xml_io_base.h"
#include "xml_io_stream.h"

template <typename... Ts>
struct xml_io_stream_name<std::tuple<Ts...>> {
  static constexpr std::string_view name = "Tuple"sv;
};

template <arts_xml_ioable T, arts_xml_ioable... Ts>
struct xml_io_stream<std::tuple<T, Ts...>> {
  constexpr static std::string_view type_name =
      xml_io_stream_name_v<std::tuple<T, Ts...>>;

  constexpr static bool all_binary =
      xml_io_binary<T> and (xml_io_binary<Ts> and ...);
  constexpr static bool all_parse =
      xml_io_parseable<T> and (xml_io_parseable<Ts> and ...);
  constexpr static bool all_same = (std::same_as<T, Ts> and ...);
  constexpr static Size N        = 1 + sizeof...(Ts);

  static void put(std::span<const std::tuple<T, Ts...>> x, bofstream* pbofs)
    requires(all_binary)
  {
    if constexpr (all_same) {
      xml_io_stream<T>::put(
          {reinterpret_cast<const T*>(x.data()), x.size() * N}, pbofs);
    } else {
      for (auto& v : x) {
        std::apply(
            [&pbofs](auto& t0, auto&... t) {
              xml_io_stream<T>::put(t0, pbofs);
              (xml_io_stream<Ts>::put(t, pbofs), ...);
            },
            v);
      }
    }
  }

  static void write(std::ostream& os,
                    const std::tuple<T, Ts...>& n,
                    bofstream* pbofs = nullptr,
                    std::string_view = ""sv) {
    if (pbofs) {
      if constexpr (all_binary)
        put({&n, 1}, pbofs);
      else
        std::apply(
            [&os, &pbofs](const T& t0, const Ts&... t) {
              xml_io_stream<T>::write(os, t0, pbofs);
              (xml_io_stream<Ts>::write(os, t, pbofs), ...);
            },
            n);
    } else {
      if constexpr (all_parse)
        std::println(os, "{:IO}", n);
      else
        std::apply(
            [&os, &pbofs](const T& t0, const Ts&... t) {
              xml_io_stream<T>::write(os, t0, pbofs);
              (xml_io_stream<Ts>::write(os, t, pbofs), ...);
            },
            n);
    }
  }

  static void get(std::span<std::tuple<T, Ts...>> x, bifstream* pbifs)
    requires(all_binary)
  {
    if constexpr (all_same) {
      xml_io_stream<T>::get({reinterpret_cast<T*>(x.data()), x.size() * N},
                            pbifs);
    } else {
      for (auto& v : x) {
        std::apply(
            [&pbifs](auto& t0, auto&... t) {
              xml_io_stream<T>::get(t0, pbifs);
              (xml_io_stream<Ts>::get(t, pbifs), ...);
            },
            v);
      }
    }
  }

  static void parse(std::span<std::tuple<T, Ts...>> x, std::istream& is)
    requires(all_parse)
  {
    if constexpr (all_same) {
      xml_io_stream<T>::parse({reinterpret_cast<T*>(x.data()), x.size() * N},
                              is);
    } else {
      for (auto& v : x) {
        std::apply(
            [&is](auto& t0, auto&... t) {
              xml_io_stream<T>::parse(t0, is);
              (xml_io_stream<Ts>::parse(t, is), ...);
            },
            v);
      }
    }
  }

  static void read(std::istream& is,
                   std::tuple<T, Ts...>& n,
                   bifstream* pbifs = nullptr) {
    if (pbifs) {
      if constexpr (all_binary)
        get({&n, 1}, pbifs);
      else
        std::apply(
            [&](T& t0, Ts&... t) {
              xml_io_stream<T>::read(is, t0, pbifs);
              (xml_io_stream<Ts>::read(is, t, pbifs), ...);
            },
            n);
    } else {
      if constexpr (all_parse)
        parse({&n, 1}, is);
      else
        std::apply(
            [&](T& t0, Ts&... t) {
              xml_io_stream<T>::read(is, t0, pbifs);
              (xml_io_stream<Ts>::read(is, t, pbifs), ...);
            },
            n);
    }
  }
};

template <typename A, typename B>
struct xml_io_stream_name<std::pair<A, B>> {
  static constexpr std::string_view name =
      xml_io_stream_name_v<std::tuple<A, B>>;
};

template <arts_xml_ioable A, arts_xml_ioable B>
struct xml_io_stream<std::pair<A, B>> {
  constexpr static std::string_view type_name =
      xml_io_stream_name_v<std::pair<A, B>>;

  constexpr static bool all_binary = xml_io_binary<A> and xml_io_binary<B>;
  constexpr static bool all_parse = xml_io_parseable<A> and xml_io_parseable<B>;
  constexpr static bool all_same  = std::same_as<A, B>;
  constexpr static Size N         = 2;

  static void put(std::span<const std::pair<A, B>> x, bofstream* pbofs)
    requires(all_binary)
  {
    if constexpr (all_same) {
      xml_io_stream<A>::put(
          {reinterpret_cast<const A*>(x.data()), x.size() * N}, pbofs);
    } else {
      for (auto& v : x) {
        xml_io_stream<A>::put(v.first, pbofs);
        xml_io_stream<B>::put(v.second, pbofs);
      }
    }
  }

  static void write(std::ostream& os,
                    const std::pair<A, B>& n,
                    bofstream* pbofs = nullptr,
                    std::string_view = ""sv) {
    if (pbofs) {
      if constexpr (all_binary)
        put({&n, 1}, pbofs);
      else {
        xml_io_stream<A>::write(os, n.first, pbofs);
        xml_io_stream<B>::write(os, n.second, pbofs);
      }
    } else {
      if constexpr (all_parse) {
        std::println(os, "{:IO}", n);
      } else {
        xml_io_stream<A>::write(os, n.first, pbofs);
        xml_io_stream<B>::write(os, n.second, pbofs);
      }
    }
  }

  static void get(std::span<std::pair<A, B>> x, bifstream* pbifs)
    requires(all_binary)
  {
    if constexpr (all_same) {
      xml_io_stream<A>::get({reinterpret_cast<A*>(x.data()), x.size() * N},
                            pbifs);
    } else {
      for (auto& v : x) {
        xml_io_stream<A>::get(v.first, pbifs);
        xml_io_stream<B>::get(v.second, pbifs);
      }
    }
  }

  static void parse(std::span<std::pair<A, B>> x, std::istream& is)
    requires(all_parse)
  {
    if constexpr (all_same) {
      xml_io_stream<A>::parse({reinterpret_cast<A*>(x.data()), x.size() * N},
                              is);
    } else {
      for (auto& v : x) {
        xml_io_stream<A>::parse(v.first, is);
        xml_io_stream<B>::parse(v.second, is);
      }
    }
  }

  static void read(std::istream& is,
                   std::pair<A, B>& n,
                   bifstream* pbifs = nullptr) {
    if (pbifs) {
      if constexpr (all_binary)
        get({&n, 1}, pbifs);
      else {
        xml_io_stream<A>::read(is, n.first, pbifs);
        xml_io_stream<B>::read(is, n.second, pbifs);
      }
    } else {
      if constexpr (all_parse)
        parse({&n, 1}, is);
      else {
        xml_io_stream<A>::read(is, n.first, pbifs);
        xml_io_stream<B>::read(is, n.second, pbifs);
      }
    }
  }
};
