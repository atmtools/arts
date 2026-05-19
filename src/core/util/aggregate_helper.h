#pragma once

#include <tuple>
#include <type_traits>

struct aggregate_init {
  template <typename T>
  constexpr operator T() const noexcept;
};

template <typename T>
concept aggregate_0 = std::is_aggregate_v<T> and requires { T{}; };

template <typename T>
concept aggregate_1 = aggregate_0<T> and requires { T{aggregate_init{}}; };

template <typename T>
concept aggregate_2 = aggregate_1<T> and requires {
  T{aggregate_init{}, aggregate_init{}};
};

template <typename T>
concept aggregate_3 = aggregate_2<T> and requires {
  T{aggregate_init{}, aggregate_init{}, aggregate_init{}};
};

template <typename T>
concept aggregate_4 = aggregate_3<T> and requires {
  T{aggregate_init{}, aggregate_init{}, aggregate_init{}, aggregate_init{}};
};

template <typename T>
concept aggregate_5 = aggregate_4<T> and requires {
  T{aggregate_init{}, aggregate_init{}, aggregate_init{}, aggregate_init{}, aggregate_init{}};
};

template <typename T>
concept aggregate_6 =
    aggregate_5<T> and requires {
      T{aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{}};
    };

template <typename T>
concept aggregate_7 =
    aggregate_6<T> and requires {
      T{aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{}};
    };

template <typename T>
concept aggregate_8 =
    aggregate_7<T> and requires {
      T{aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{}};
    };

template <typename T>
concept aggregate_9 =
    aggregate_8<T> and requires {
      T{aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{}};
    };

template <typename T>
concept aggregate_10 =
    aggregate_9<T> and requires {
      T{aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{},
        aggregate_init{}};
    };

template <typename T>
concept aggregate_11 = aggregate_10<T> and requires {
  T{aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{}};
};

template <typename T>
concept aggregate_12 = aggregate_11<T> and requires {
  T{aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{}};
};

template <typename T>
concept aggregate_13 = aggregate_12<T> and requires {
  T{aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{}};
};

template <typename T>
concept aggregate_14 = aggregate_13<T> and requires {
  T{aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{}};
};

template <typename T>
concept aggregate_15 = aggregate_14<T> and requires {
  T{aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{}};
};

template <typename T>
concept aggregate_16 = aggregate_15<T> and requires {
  T{aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{}};
};

template <typename T>
concept aggregate_17 = aggregate_16<T> and requires {
  T{aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{}};
};

template <typename T>
concept aggregate_18 = aggregate_17<T> and requires {
  T{aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{}};
};

template <typename T>
concept aggregate_19 = aggregate_18<T> and requires {
  T{aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{}};
};

template <typename T>
concept aggregate_20 = aggregate_19<T> and requires {
  T{aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{},
    aggregate_init{}};
};

template <typename T>
concept aggregate_0_exact = aggregate_0<T> and not aggregate_1<T>;

template <typename T>
concept aggregate_1_exact = aggregate_1<T> and not aggregate_2<T>;

template <typename T>
concept aggregate_2_exact = aggregate_2<T> and not aggregate_3<T>;

template <typename T>
concept aggregate_3_exact = aggregate_3<T> and not aggregate_4<T>;

template <typename T>
concept aggregate_4_exact = aggregate_4<T> and not aggregate_5<T>;

template <typename T>
concept aggregate_5_exact = aggregate_5<T> and not aggregate_6<T>;

template <typename T>
concept aggregate_6_exact = aggregate_6<T> and not aggregate_7<T>;

template <typename T>
concept aggregate_7_exact = aggregate_7<T> and not aggregate_8<T>;

template <typename T>
concept aggregate_8_exact = aggregate_8<T> and not aggregate_9<T>;

template <typename T>
concept aggregate_9_exact = aggregate_9<T> and not aggregate_10<T>;

template <typename T>
concept aggregate_10_exact = aggregate_10<T> and not aggregate_11<T>;

template <typename T>
concept aggregate_11_exact = aggregate_11<T> and not aggregate_12<T>;

template <typename T>
concept aggregate_12_exact = aggregate_12<T> and not aggregate_13<T>;

template <typename T>
concept aggregate_13_exact = aggregate_13<T> and not aggregate_14<T>;

template <typename T>
concept aggregate_14_exact = aggregate_14<T> and not aggregate_15<T>;

template <typename T>
concept aggregate_15_exact = aggregate_15<T> and not aggregate_16<T>;

template <typename T>
concept aggregate_16_exact = aggregate_16<T> and not aggregate_17<T>;

template <typename T>
concept aggregate_17_exact = aggregate_17<T> and not aggregate_18<T>;

template <typename T>
concept aggregate_18_exact = aggregate_18<T> and not aggregate_19<T>;

template <typename T>
concept aggregate_19_exact = aggregate_19<T> and not aggregate_20<T>;

template <typename T>
concept aggregate_20_exact = aggregate_20<T>;

auto as_tuple(aggregate_0_exact auto&) { return std::tuple<>(); }

auto as_tuple(aggregate_1_exact auto& x) {
  auto&& [a] = x;
  return std::tie(a);
}

auto as_tuple(aggregate_2_exact auto& x) {
  auto&& [a, b] = x;
  return std::tie(a, b);
}

auto as_tuple(aggregate_3_exact auto& x) {
  auto&& [a, b, c] = x;
  return std::tie(a, b, c);
}

auto as_tuple(aggregate_4_exact auto& x) {
  auto&& [a, b, c, d] = x;
  return std::tie(a, b, c, d);
}

auto as_tuple(aggregate_5_exact auto& x) {
  auto&& [a, b, c, d, e] = x;
  return std::tie(a, b, c, d, e);
}

auto as_tuple(aggregate_6_exact auto& x) {
  auto&& [a, b, c, d, e, f] = x;
  return std::tie(a, b, c, d, e, f);
}

auto as_tuple(aggregate_7_exact auto& x) {
  auto&& [a, b, c, d, e, f, g] = x;
  return std::tie(a, b, c, d, e, f, g);
}

auto as_tuple(aggregate_8_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h] = x;
  return std::tie(a, b, c, d, e, f, g, h);
}

auto as_tuple(aggregate_9_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i] = x;
  return std::tie(a, b, c, d, e, f, g, h, i);
}

auto as_tuple(aggregate_10_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j);
}

auto as_tuple(aggregate_11_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k);
}

auto as_tuple(aggregate_12_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l);
}

auto as_tuple(aggregate_13_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m);
}

auto as_tuple(aggregate_14_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n);
}

auto as_tuple(aggregate_15_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o);
}

auto as_tuple(aggregate_16_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p);
}

auto as_tuple(aggregate_17_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q);
}

auto as_tuple(aggregate_18_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r);
}

auto as_tuple(aggregate_19_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s);
}

auto as_tuple(aggregate_20_exact auto& x) {
  auto&& [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t] = x;
  return std::tie(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t);
}

template <typename T>
concept arts_aggregate = requires(T a) { as_tuple(a); };
