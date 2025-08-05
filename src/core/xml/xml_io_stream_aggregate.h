#pragma once

#include <concepts>
#include <type_traits>

#include "xml_io_base.h"
#include "xml_io_stream.h"

template <typename T>
concept aggregate_0 = std::is_aggregate_v<T> and requires { T{}; };

template <typename T>
concept aggregate_1 = aggregate_0<T> and requires { T({}); };

template <typename T>
concept aggregate_2 = aggregate_1<T> and requires { T({}, {}); };

template <typename T>
concept aggregate_3 = aggregate_2<T> and requires { T({}, {}, {}); };

template <typename T>
concept aggregate_4 = aggregate_3<T> and requires { T({}, {}, {}, {}); };

template <typename T>
concept aggregate_5 = aggregate_4<T> and requires { T({}, {}, {}, {}, {}); };

template <typename T>
concept aggregate_6 =
    aggregate_5<T> and requires { T({}, {}, {}, {}, {}, {}); };

template <typename T>
concept aggregate_7 =
    aggregate_6<T> and requires { T({}, {}, {}, {}, {}, {}, {}); };

template <typename T>
concept aggregate_8 =
    aggregate_7<T> and requires { T({}, {}, {}, {}, {}, {}, {}, {}); };

template <typename T>
concept aggregate_9 =
    aggregate_8<T> and requires { T({}, {}, {}, {}, {}, {}, {}, {}, {}); };

template <typename T>
concept aggregate_10 =
    aggregate_9<T> and requires { T({}, {}, {}, {}, {}, {}, {}, {}, {}, {}); };

template <typename T>
concept aggregate_11 = aggregate_10<T> and requires {
  T({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {});
};

template <typename T>
concept aggregate_12 = aggregate_11<T> and requires {
  T({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {});
};

template <typename T>
concept aggregate_13 = aggregate_12<T> and requires {
  T({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {});
};

template <typename T>
concept aggregate_14 = aggregate_13<T> and requires {
  T({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {});
};

template <typename T>
concept aggregate_15 = aggregate_14<T> and requires {
  T({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {});
};

template <typename T>
concept aggregate_16 = aggregate_15<T> and requires {
  T({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {});
};

template <typename T>
concept aggregate_17 = aggregate_16<T> and requires {
  T({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {});
};

template <typename T>
concept aggregate_18 = aggregate_17<T> and requires {
  T({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {});
};

template <typename T>
concept aggregate_19 = aggregate_18<T> and requires {
  T({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {});
};

template <typename T>
concept aggregate_20 = aggregate_19<T> and requires {
  T({},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {},
    {});
};

auto as_tuple(aggregate_0 auto&) { return std::tuple<>(); }
auto as_tuple(aggregate_1 auto& x) {
  auto&& [a] = x;
  return std::tie(a);
}

auto as_tuple(aggregate_2 auto& x) {
  auto&& [a, b] = x;
  return std::tie(a, b);
}

auto as_tuple(aggregate_3 auto& x) {
  auto&& [a, b, c] = x;
  return std::tie(a, b, c);
}

auto as_tuple(aggregate_4 auto& x) {
  auto&& [a, b, c, d] = x;
  return std::tie(a, b, c, d);
}

auto as_tuple(aggregate_5 auto& x) {
  auto&& [a, b, c, d, e] = x;
  return std::tie(a, b, c, d, e);
}

auto as_tuple(aggregate_6 auto& x) {
  auto&& [a, b, c, d, e, f] = x;
  return std::tie(a, b, c, d, e, f);
}

auto as_tuple(aggregate_7 auto& x) {
  auto&& [a, b, c, d, e, f, g] = x;
  return std::tie(a, b, c, d, e, f, g);
}

auto as_tuple(aggregate_8 auto& x) {
  auto&& [a, b, c, d, e, f, g, h] = x;
  return std::tie(a, b, c, d, e, f, g, h);
}

auto as_tuple(aggregate_9 auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i] = x;
  return std::tie(a, b, c, d, e, f, g, h, i);
}

auto as_tuple(aggregate_10 auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j);
}

auto as_tuple(aggregate_11 auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k);
}

auto as_tuple(aggregate_12 auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l);
}

auto as_tuple(aggregate_13 auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m);
}

auto as_tuple(aggregate_14 auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n);
}

auto as_tuple(aggregate_15 auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o);
}

auto as_tuple(aggregate_16 auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p);
}

auto as_tuple(aggregate_17 auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q);
}

auto as_tuple(aggregate_18 auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r);
}

auto as_tuple(aggregate_19 auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s);
}

auto as_tuple(aggregate_20 auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t);
}

template <typename T>
struct xml_io_stream_aggregate {
  constexpr static bool value = false;
};

template <typename T>
constexpr bool xml_io_stream_aggregate_v = xml_io_stream_aggregate<T>::value;

template <typename T>
concept xml_io_aggregratable =
    xml_io_stream_aggregate_v<T> and requires(T a) { as_tuple(a); };

template <xml_io_aggregratable T>
struct xml_io_stream<T> {
  constexpr static std::string_view type_name = xml_io_stream_name_v<T>;

  using _lambda = decltype([](auto& v) { return as_tuple(v); });
  using ctup_t  = std::invoke_result_t<_lambda, const T&>;
  using mtup_t  = std::invoke_result_t<_lambda, T&>;
  using inner   = xml_io_stream<mtup_t>;

  static void write(std::ostream& os,
                    const T& t,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

    std::apply(
        [&os, &pbofs]<typename... Ts>(const Ts&... v) {
          (xml_io_stream<Ts>::write(os, v, pbofs), ...);
        },
        as_tuple(t));

  tag.write_to_end_stream(os);
  }

  static void read(std::istream& is, T& t, bifstream* pbifs = nullptr) try {
    XMLTag tag{};
    tag.read_from_stream(is);
    tag.check_name(type_name);

    std::apply([&is, &pbifs]<typename... Ts>(
                   Ts&... v) { (xml_io_stream<Ts>::read(is, v, pbifs), ...); },
               as_tuple(t));

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::format("Error reading {}:\n{}", type_name, e.what()));
  }
};
