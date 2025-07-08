#pragma once

#include <format_tags.h>

#include <concepts>
#include <tuple>

#include "xml_io_base.h"
#include "xml_io_stream.h"

template <typename... Ts>
struct xml_io_stream_name<std::tuple<Ts...>> {
  static constexpr std::string_view name = "Tuple"sv;
};

template <typename... Ts>
  requires((arts_xml_ioable<std::remove_cvref_t<Ts>> and ...))
struct xml_io_stream<std::tuple<Ts...>> {
  template <typename U>
  using Decayed = std::remove_cvref_t<U>;

  constexpr static std::string_view type_name =
      xml_io_stream_name_v<std::tuple<Ts...>>;

  constexpr static bool all_binary = (xml_io_binary<Decayed<Ts>> and ...);
  constexpr static bool all_parse  = (xml_io_parseable<Decayed<Ts>> and ...);

  static void put(std::span<const std::tuple<Ts...>> x, bofstream* pbofs)
    requires(all_binary)
  {
    for (auto& v : x) {
      std::apply(
          [&pbofs](auto&... args) {
            (xml_io_stream<Decayed<decltype(args)>>::put({&args, 1}, pbofs),
             ...);
          },
          v);
    }
  }

  static void write(std::ostream& os,
                    const std::tuple<Ts...>& n,
                    bofstream* pbofs = nullptr,
                    std::string_view = ""sv) {
    if (pbofs) {
      if constexpr (all_binary)
        put({&n, 1}, pbofs);
      else
        std::apply(
            [&os, &pbofs](const auto&... args) {
              (xml_io_stream<Decayed<decltype(args)>>::write(os, args, pbofs),
               ...);
            },
            n);
    } else {
      if constexpr (all_parse)
        std::println(os, "{:IO}", n);
      else
        std::apply(
            [&os, &pbofs](const auto&... args) {
              (xml_io_stream<Decayed<decltype(args)>>::write(os, args, pbofs),
               ...);
            },
            n);
    }
  }

  static void get(std::span<std::tuple<Ts...>> x, bifstream* pbifs)
    requires(all_binary)
  {
    for (auto& v : x) {
      std::apply(
          [&pbifs](auto&... args) {
            (xml_io_stream<Decayed<decltype(args)>>::get(args, pbifs), ...);
          },
          v);
    }
  }

  static void parse(std::span<std::tuple<Ts...>> x, std::istream& is)
    requires(all_parse)
  {
    for (auto& v : x) {
      std::apply(
          [&is](auto&... args) {
            (xml_io_stream<Decayed<decltype(args)>>::parse(args, is), ...);
          },
          v);
    }
  }

  static void read(std::istream& is,
                   std::tuple<Ts...>& n,
                   bifstream* pbifs = nullptr) try {
    if (pbifs) {
      if constexpr (all_binary)
        get({&n, 1}, pbifs);
      else
        std::apply(
            [&](auto&... args) {
              (xml_io_stream<Decayed<decltype(args)>>::read(is, args, pbifs),
               ...);
            },
            n);
    } else {
      if constexpr (all_parse)
        parse({&n, 1}, is);
      else
        std::apply(
            [&](auto&... args) {
              (xml_io_stream<Decayed<decltype(args)>>::read(is, args, pbifs),
               ...);
            },
            n);
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::format("Error reading {}<{:,}>:\n{}",
                    type_name,
                    std::array{xml_io_stream<Ts>::type_name...},
                    e.what()));
  }
};

template <typename A, typename B>
struct xml_io_stream_name<std::pair<A, B>> {
  static constexpr std::string_view name =
      xml_io_stream_name_v<std::tuple<A, B>>;
};

template <typename A, typename B>
  requires(arts_xml_ioable<std::remove_cvref_t<A>> and
           arts_xml_ioable<std::remove_cvref_t<B>>)
struct xml_io_stream<std::pair<A, B>> {
  template <typename U>
  using Decayed = std::remove_cvref_t<U>;

  constexpr static std::string_view type_name =
      xml_io_stream_name_v<std::pair<A, B>>;

  constexpr static bool all_binary =
      xml_io_binary<Decayed<A>> and xml_io_binary<Decayed<B>>;
  constexpr static bool all_parse =
      xml_io_parseable<Decayed<A>> and xml_io_parseable<Decayed<B>>;

  static void put(std::span<const std::pair<A, B>> x, bofstream* pbofs)
    requires(all_binary)
  {
    for (auto& v : x) {
      xml_io_stream<Decayed<A>>::put({&v.first, 1}, pbofs);
      xml_io_stream<Decayed<B>>::put({&v.second, 1}, pbofs);
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
        xml_io_stream<Decayed<A>>::write(os, n.first, pbofs);
        xml_io_stream<Decayed<B>>::write(os, n.second, pbofs);
      }
    } else {
      if constexpr (all_parse) {
        std::println(os, "{:IO}", n);
      } else {
        xml_io_stream<Decayed<A>>::write(os, n.first, pbofs);
        xml_io_stream<Decayed<B>>::write(os, n.second, pbofs);
      }
    }
  }

  static void get(std::span<std::pair<A, B>> x, bifstream* pbifs)
    requires(all_binary)
  {
    for (auto& v : x) {
      xml_io_stream<Decayed<A>>::get({&v.first, 1}, pbifs);
      xml_io_stream<Decayed<B>>::get({&v.second, 1}, pbifs);
    }
  }

  static void parse(std::span<std::pair<A, B>> x, std::istream& is)
    requires(all_parse)
  {
    for (auto& v : x) {
      xml_io_stream<Decayed<A>>::parse({&v.first, 1}, is);
      xml_io_stream<Decayed<B>>::parse({&v.second, 1}, is);
    }
  }

  static void read(std::istream& is,
                   std::pair<A, B>& n,
                   bifstream* pbifs = nullptr) {
    if (pbifs) {
      if constexpr (all_binary)
        get({&n, 1}, pbifs);
      else {
        xml_io_stream<Decayed<A>>::read(is, n.first, pbifs);
        xml_io_stream<Decayed<B>>::read(is, n.second, pbifs);
      }
    } else {
      if constexpr (all_parse)
        parse({&n, 1}, is);
      else {
        xml_io_stream<Decayed<A>>::read(is, n.first, pbifs);
        xml_io_stream<Decayed<B>>::read(is, n.second, pbifs);
      }
    }
  }
};
