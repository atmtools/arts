#include <configtypes.h>
#include <lagrange_interp.h>
#include <matpack.h>
#include <matpack_mdspan_helpers_vector.h>
#include <rng.h>

#include <cstdlib>

#include "matpack_mdspan_elemwise_mditer.h"
#include "time_test_util.h"

namespace {
template <lagrange_interp::transformer t>
consteval std::string_view transform_name() {
  if constexpr (std::same_as<t, lagrange_interp::identity>) {
    return "identity"sv;
  } else if constexpr (std::same_as<t, lagrange_interp::loncross>) {
    return "loncross"sv;
  }

  return "<unknown>"sv;
}

bool close(const Numeric& a, const Numeric& b) {
  return nonstd::abs(a - b) < 1e-6;
}

bool all_close(const std::vector<Numeric>& a) {
  if (a.empty()) return true;

  const Numeric ref = a.front();
  for (const auto& v : a) {
    if (not close(v, ref)) return false;
  }

  return true;
}

template <Size N,
          lagrange_interp::transformer Transformer1,
          lagrange_interp::transformer Transformer2>
void test_fixed() {
  using namespace lagrange_interp;

  constexpr auto name1 = transform_name<Transformer1>();
  constexpr auto name2 = transform_name<Transformer2>();

  const Vector grid1{-120, -90, -60, -30, 0, 30, 60, 90, 120};
  const Vector grid2{120, 90, 60, 30, 0, -30, -60, -90, -120};

  const Vector coords1{nlinspace(-180, 180, 300)};
  const Vector coords2{nlinspace(180, -180, 300)};

  const Matrix data = nlinspace(1, 2, grid1.size() * grid2.size())
                          .reshape(grid1.size(), grid2.size());

  std::vector<Numeric> results_interp;
  results_interp.reserve(4);

  std::vector<Numeric> results_reinterp;
  results_reinterp.reserve(4);

  std::vector<Numeric> results_flat_interp;
  results_flat_interp.reserve(4);

  std::vector<lag_t<N, Transformer1>> lags11(coords1.size());
  std::vector<lag_t<N, Transformer1>> lags12(coords1.size());
  std::vector<lag_t<N, Transformer2>> lags21(coords2.size());
  std::vector<lag_t<N, Transformer2>> lags22(coords2.size());

  {
    test_timer_t timer(std::format(
        "make_lags fixed 1-by-1 linear order {}, {}-{}", N, name1, name2));

    for (Size i = 0; i < coords1.size(); ++i) {
      lags11[i] = lag_t<N, Transformer1>(grid1, coords1[i], ascending_grid_t{});
      lags21[i] =
          lag_t<N, Transformer2>(grid2, coords2[i], descending_grid_t{});
    }
  }

  {
    test_timer_t timer(std::format(
        "make_lags fixed all-in-one linear order {}, {}-{}", N, name1, name2));

    lags12 = make_lags<N, Transformer1>(grid1, coords1, 0.0, "grid1");
    lags22 = make_lags<N, Transformer2>(grid2, coords2, 0.0, "grid2");
  }

  // Test direct interp 1
  {
    Numeric& res = results_interp.emplace_back();

    test_timer_t timer(std::format(
        "interp fixed 1-by-1 linear order {}, {}-{}", N, name1, name2));

    for (Size i = 0; i < coords1.size(); ++i) {
      for (Size j = 0; j < coords2.size(); ++j) {
        res += interp(data, lags11[i], lags21[j]);
      }
    }
  }

  // Test direct interp 2
  {
    Numeric& res = results_interp.emplace_back();

    test_timer_t timer(std::format(
        "interp fixed all-in-one linear order {}, {}-{}", N, name1, name2));

    for (Size i = 0; i < coords1.size(); ++i) {
      for (Size j = 0; j < coords2.size(); ++j) {
        res += interp(data, lags12[i], lags22[j]);
      }
    }
  }

  // Test direct reinterp 1
  {
    Numeric& res = results_reinterp.emplace_back();

    test_timer_t timer(std::format(
        "reinterp fixed 1-by-1 linear order {}, {}-{}", N, name1, name2));

    res = sum(reinterp(data, lags11, lags21));
  }

  // Test direct reinterp 2
  {
    Numeric& res = results_reinterp.emplace_back();

    test_timer_t timer(std::format(
        "reinterp fixed all-in-one linear order {}, {}-{}", N, name1, name2));

    res = sum(reinterp(data, lags12, lags22));
  }

  // Test direct flat_interp 1
  {
    Numeric& res = results_flat_interp.emplace_back();

    test_timer_t timer(std::format(
        "flat_interp fixed 1-by-1 linear order {}, {}-{}", N, name1, name2));

    res = sum(flat_interp(data, lags11, lags21));
  }

  // Test direct flat_interp 2
  {
    Numeric& res = results_flat_interp.emplace_back();

    test_timer_t timer(
        std::format("flat_interp fixed all-in-one linear order {}, {}-{}",
                    N,
                    name1,
                    name2));

    res = sum(flat_interp(data, lags12, lags22));
  }

  // Test indirect interp 1
  {
    Numeric& res = results_interp.emplace_back();

    test_timer_t timer(std::format(
        "interp-itw fixed 1-by-1 linear order {}, {}-{}", N, name1, name2));

    for (Size i = 0; i < coords1.size(); ++i) {
      for (Size j = 0; j < coords2.size(); ++j) {
        res += interp(
            data, interpweights(lags11[i], lags21[j]), lags11[i], lags21[j]);
      }
    }
  }

  // Test indirect interp 2
  {
    Numeric& res = results_interp.emplace_back();

    test_timer_t timer(std::format(
        "interp-itw fixed all-in-one linear order {}, {}-{}", N, name1, name2));

    for (Size i = 0; i < coords1.size(); ++i) {
      for (Size j = 0; j < coords2.size(); ++j) {
        res += interp(
            data, interpweights(lags12[i], lags22[j]), lags12[i], lags22[j]);
      }
    }
  }

  // Test indirect reinterp 1
  {
    Numeric& res = results_reinterp.emplace_back();

    test_timer_t timer(std::format(
        "reinterp-itw fixed 1-by-1 linear order {}, {}-{}", N, name1, name2));

    res = sum(reinterp(data, reinterpweights(lags11, lags21), lags11, lags21));
  }

  // Test indirect reinterp 2
  {
    Numeric& res = results_reinterp.emplace_back();

    test_timer_t timer(
        std::format("reinterp-itw fixed all-in-one linear order {}, {}-{}",
                    N,
                    name1,
                    name2));

    res = sum(reinterp(data, reinterpweights(lags12, lags22), lags12, lags22));
  }

  // Test indirect flat_interp 1
  {
    Numeric& res = results_flat_interp.emplace_back();

    test_timer_t timer(
        std::format("flat_interp-itw fixed 1-by-1 linear order {}, {}-{}",
                    N,
                    name1,
                    name2));

    res = sum(
        flat_interp(data, flat_interpweights(lags11, lags21), lags11, lags21));
  }

  // Test indirect flat_interp 2
  {
    Numeric& res = results_flat_interp.emplace_back();

    test_timer_t timer(
        std::format("flat_interp-itw fixed all-in-one linear order {}, {}-{}",
                    N,
                    name1,
                    name2));

    res = sum(
        flat_interp(data, flat_interpweights(lags12, lags22), lags12, lags22));
  }

  if (not all_close(results_interp) or not all_close(results_reinterp) or
      not all_close(results_flat_interp)) {
    std::println(R"(Results are not close for {}-{} with N={}

interp:      {:B,}
reinterp:    {:B,},
flat_interp: {:B,})",
                 name1,
                 name2,
                 N,
                 results_interp,
                 results_reinterp,
                 results_flat_interp);
  }
}

template <lagrange_interp::transformer Transformer1,
          lagrange_interp::transformer Transformer2>
void test_runtime(Size N) {
  using namespace lagrange_interp;

  constexpr auto name1 = transform_name<Transformer1>();
  constexpr auto name2 = transform_name<Transformer2>();

  const Vector grid1{-120, -90, -60, -30, 0, 30, 60, 90, 120};
  const Vector grid2{120, 90, 60, 30, 0, -30, -60, -90, -120};

  const Vector coords1{nlinspace(-180, 180, 300)};
  const Vector coords2{nlinspace(180, -180, 300)};

  const Matrix data = nlinspace(1, 2, grid1.size() * grid2.size())
                          .reshape(grid1.size(), grid2.size());

  std::vector<Numeric> results_interp;
  results_interp.reserve(4);

  std::vector<Numeric> results_reinterp;
  results_reinterp.reserve(4);

  std::vector<Numeric> results_flat_interp;
  results_flat_interp.reserve(4);

  std::vector<lag_t<-1, Transformer1>> lags11(coords1.size());
  std::vector<lag_t<-1, Transformer1>> lags12(coords1.size());
  std::vector<lag_t<-1, Transformer2>> lags21(coords2.size());
  std::vector<lag_t<-1, Transformer2>> lags22(coords2.size());

  {
    test_timer_t timer(std::format(
        "make_lags dynamic 1-by-1 linear order {}, {}-{}", N, name1, name2));

    for (Size i = 0; i < coords1.size(); ++i) {
      lags11[i] =
          lag_t<-1, Transformer1>(grid1, coords1[i], N, ascending_grid_t{});
      lags21[i] =
          lag_t<-1, Transformer2>(grid2, coords2[i], N, descending_grid_t{});
    }
  }

  {
    test_timer_t timer(
        std::format("make_lags dynamic all-in-one linear order {}, {}-{}",
                    N,
                    name1,
                    name2));

    lags12 = make_lags<Transformer1>(grid1, coords1, N, 0.0, "grid1");
    lags22 = make_lags<Transformer2>(grid2, coords2, N, 0.0, "grid2");
  }

  // Test direct interp 1
  {
    Numeric& res = results_interp.emplace_back();

    test_timer_t timer(std::format(
        "interp dynamic 1-by-1 linear order {}, {}-{}", N, name1, name2));

    for (Size i = 0; i < coords1.size(); ++i) {
      for (Size j = 0; j < coords2.size(); ++j) {
        res += interp(data, lags11[i], lags21[j]);
      }
    }
  }

  // Test direct interp 2
  {
    Numeric& res = results_interp.emplace_back();

    test_timer_t timer(std::format(
        "interp dynamic all-in-one linear order {}, {}-{}", N, name1, name2));

    for (Size i = 0; i < coords1.size(); ++i) {
      for (Size j = 0; j < coords2.size(); ++j) {
        res += interp(data, lags12[i], lags22[j]);
      }
    }
  }

  // Test direct reinterp 1
  {
    Numeric& res = results_reinterp.emplace_back();

    test_timer_t timer(std::format(
        "reinterp dynamic 1-by-1 linear order {}, {}-{}", N, name1, name2));

    res = sum(reinterp(data, lags11, lags21));
  }

  // Test direct reinterp 2
  {
    Numeric& res = results_reinterp.emplace_back();

    test_timer_t timer(std::format(
        "reinterp dynamic all-in-one linear order {}, {}-{}", N, name1, name2));

    res = sum(reinterp(data, lags12, lags22));
  }

  // Test direct flat_interp 1
  {
    Numeric& res = results_flat_interp.emplace_back();

    test_timer_t timer(std::format(
        "flat_interp dynamic 1-by-1 linear order {}, {}-{}", N, name1, name2));

    res = sum(flat_interp(data, lags11, lags21));
  }

  // Test direct flat_interp 2
  {
    Numeric& res = results_flat_interp.emplace_back();

    test_timer_t timer(
        std::format("flat_interp dynamic all-in-one linear order {}, {}-{}",
                    N,
                    name1,
                    name2));

    res = sum(flat_interp(data, lags12, lags22));
  }

  // Test indirect interp 1
  {
    Numeric& res = results_interp.emplace_back();

    test_timer_t timer(std::format(
        "interp-itw dynamic 1-by-1 linear order {}, {}-{}", N, name1, name2));

    for (Size i = 0; i < coords1.size(); ++i) {
      for (Size j = 0; j < coords2.size(); ++j) {
        res += interp(
            data, interpweights(lags11[i], lags21[j]), lags11[i], lags21[j]);
      }
    }
  }

  // Test indirect interp 2
  {
    Numeric& res = results_interp.emplace_back();

    test_timer_t timer(
        std::format("interp-itw dynamic all-in-one linear order {}, {}-{}",
                    N,
                    name1,
                    name2));

    for (Size i = 0; i < coords1.size(); ++i) {
      for (Size j = 0; j < coords2.size(); ++j) {
        res += interp(
            data, interpweights(lags12[i], lags22[j]), lags12[i], lags22[j]);
      }
    }
  }

  // Test indirect reinterp 1
  {
    Numeric& res = results_reinterp.emplace_back();

    test_timer_t timer(std::format(
        "reinterp-itw dynamic 1-by-1 linear order {}, {}-{}", N, name1, name2));

    res = sum(reinterp(data, reinterpweights(lags11, lags21), lags11, lags21));
  }

  // Test indirect reinterp 2
  {
    Numeric& res = results_reinterp.emplace_back();

    test_timer_t timer(
        std::format("reinterp-itw dynamic all-in-one linear order {}, {}-{}",
                    N,
                    name1,
                    name2));

    res = sum(reinterp(data, reinterpweights(lags12, lags22), lags12, lags22));
  }

  // Test indirect flat_interp 1
  {
    Numeric& res = results_flat_interp.emplace_back();

    test_timer_t timer(
        std::format("flat_interp-itw dynamic 1-by-1 linear order {}, {}-{}",
                    N,
                    name1,
                    name2));

    res = sum(
        flat_interp(data, flat_interpweights(lags11, lags21), lags11, lags21));
  }

  // Test indirect flat_interp 2
  {
    Numeric& res = results_flat_interp.emplace_back();

    test_timer_t timer(
        std::format("flat_interp-itw dynamic all-in-one linear order {}, {}-{}",
                    N,
                    name1,
                    name2));

    res = sum(
        flat_interp(data, flat_interpweights(lags12, lags22), lags12, lags22));
  }

  if (not all_close(results_interp) or not all_close(results_reinterp) or
      not all_close(results_flat_interp)) {
    std::println(R"(Results are not close for {}-{} with N={}

interp:      {:B,}
reinterp:    {:B,},
flat_interp: {:B,})",
                 name1,
                 name2,
                 N,
                 results_interp,
                 results_reinterp,
                 results_flat_interp);
  }
}
}  // namespace

int main() {
  for (Index i = 0; i < 5; ++i) {
    test_fixed<0, lagrange_interp::identity, lagrange_interp::identity>();
    test_fixed<1, lagrange_interp::identity, lagrange_interp::identity>();
    test_fixed<2, lagrange_interp::identity, lagrange_interp::identity>();
    test_fixed<3, lagrange_interp::identity, lagrange_interp::identity>();
    test_fixed<4, lagrange_interp::identity, lagrange_interp::identity>();

    test_fixed<0, lagrange_interp::identity, lagrange_interp::loncross>();
    test_fixed<1, lagrange_interp::identity, lagrange_interp::loncross>();
    test_fixed<2, lagrange_interp::identity, lagrange_interp::loncross>();
    test_fixed<3, lagrange_interp::identity, lagrange_interp::loncross>();
    test_fixed<4, lagrange_interp::identity, lagrange_interp::loncross>();

    test_fixed<0, lagrange_interp::loncross, lagrange_interp::identity>();
    test_fixed<1, lagrange_interp::loncross, lagrange_interp::identity>();
    test_fixed<2, lagrange_interp::loncross, lagrange_interp::identity>();
    test_fixed<3, lagrange_interp::loncross, lagrange_interp::identity>();
    test_fixed<4, lagrange_interp::loncross, lagrange_interp::identity>();

    test_fixed<0, lagrange_interp::loncross, lagrange_interp::loncross>();
    test_fixed<1, lagrange_interp::loncross, lagrange_interp::loncross>();
    test_fixed<2, lagrange_interp::loncross, lagrange_interp::loncross>();
    test_fixed<3, lagrange_interp::loncross, lagrange_interp::loncross>();
    test_fixed<4, lagrange_interp::loncross, lagrange_interp::loncross>();

    for (Size N = 0; N < 5; ++N) {
      test_runtime<lagrange_interp::identity, lagrange_interp::identity>(N);
      test_runtime<lagrange_interp::identity, lagrange_interp::loncross>(N);
      test_runtime<lagrange_interp::loncross, lagrange_interp::identity>(N);
      test_runtime<lagrange_interp::loncross, lagrange_interp::loncross>(N);
    }
  }

  print_time_points();
}
