#include <configtypes.h>
#include <lagrange_interp.h>
#include <matpack.h>
#include <matpack_mdspan_helpers_vector.h>
#include <rng.h>

#include <cstdlib>

#include "time_test_util.h"

namespace {
bool allclose(Numeric x, auto... xs) {
  return ((nonstd::abs(x - xs) <
           x * 100 * std::numeric_limits<Numeric>::epsilon()) and
          ...);
}

template <Size N>
void test_fixed(bool sorted) {
  const Vector grid{1, 2, 3, 4, 5, 6};

  const Matrix data =
      matpack::uniform_grid(-180, grid.size() * grid.size(), 0.1)
          .reshape(grid.size(), grid.size());

  Numeric sumf0{}, sumf1{}, sumf2{}, sumf3{}, sumf4{};

  auto r  = random_numbers(10000, -5.0, 10.0);
  auto r2 = random_numbers(r.size(), -5.0, 10.0);

  if (sorted) {
    std::sort(r.begin(), r.end());
    std::sort(r2.begin(), r2.end());
  }

  using newT = lagrange_interp::lag_t<N>;
  std::vector<newT> flag1(r.size());
  std::vector<newT> flag2(r.size());

  {
    test_timer_t t{
        std::format("lag fixed {} {}", N, sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      flag1[i] = newT(grid, r[i]);
      flag2[i] = newT(grid, r2[i]);
    }
  }

  {
    test_timer_t t{std::format(
        "interp fixed {} {}", N, sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      sumf0 += interp(data, flag1[i], flag2[i]);
    }
  }

  {
    test_timer_t t{std::format(
        "interp-itw fixed {} {}", N, sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      const auto itw  = interpweights(flag1[i], flag2[i]);
      sumf3          += interp(data, itw, flag1[i], flag2[i]);
    }
  }

  {
    test_timer_t t{std::format(
        "flat_interp fixed {} {}", N, sorted ? "sorted"sv : "unsorted"sv)};
    sumf1 = sum(flat_interp(data, flag1, flag2));
  }

  {
    test_timer_t t{std::format("flat_interp-itw fixed {} old-N/A {}",
                               N,
                               sorted ? "sorted"sv : "unsorted"sv)};
  }

  {
    test_timer_t t{std::format(
        "flat_interp-itw fixed {} {}", N, sorted ? "sorted"sv : "unsorted"sv)};
    const auto itw = flat_interpweights(flag1, flag2);
    sumf4          = sum(flat_interp(data, flag1, flag2));
  }

  {
    test_timer_t t{std::format(
        "flat_interp fixed-vec {} {}", N, sorted ? "sorted"sv : "unsorted"sv)};
    const auto l1 = lagrange_interp::make_lags<N>(grid, r, 0.0);
    const auto l2 = lagrange_interp::make_lags<N>(grid, r2, 0.0);
    sumf2         = sum(flat_interp(data, l1, l2));
  }

  ARTS_USER_ERROR_IF(allclose(sumf0, sumf1, sumf2, sumf3, sumf4),
                     "Results of old and interpolation do not match!\n"
                     "new: {0} != {1} != {2} != {3} != {4}",
                     sumf0,
                     sumf1,
                     sumf2,
                     sumf3,
                     sumf4);
}

template <Size N>
void test_fixed_cyclic(bool sorted) {
  const Vector grid{-180.0, -100.0, -50.0, 0.0, 50.0, 100.0, 179.0};

  const Matrix data =
      matpack::uniform_grid(-180, grid.size() * grid.size(), 0.1)
          .reshape(grid.size(), grid.size());

  Numeric sumf0{}, sumf1{}, sumf2{}, sumf3{}, sumf4{};

  auto r  = random_numbers(10000, -500.0, 500.0);
  auto r2 = random_numbers(r.size(), -500.0, 500.0);

  if (sorted) {
    std::sort(r.begin(), r.end());
    std::sort(r2.begin(), r2.end());
  }

  using newT = lagrange_interp::lag_t<N, lagrange_interp::loncross>;
  std::vector<newT> flag1(r.size());
  std::vector<newT> flag2(r.size());

  {
    test_timer_t t{std::format(
        "lag fixed-cyclic {} {}", N, sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      flag1[i] = newT(grid, r[i]);
      flag2[i] = newT(grid, r2[i]);
    }
  }

  {
    test_timer_t t{std::format(
        "interp fixed-cyclic {} {}", N, sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      sumf0 += interp(data, flag1[i], flag2[i]);
    }
  }

  {
    test_timer_t t{std::format("interp-itw fixed-cyclic {} {}",
                               N,
                               sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      auto itw  = interpweights(flag1[i], flag2[i]);
      sumf3    += interp(data, itw, flag1[i], flag2[i]);
    }
  }

  {
    test_timer_t t{std::format("flat_interp fixed-cyclic {} {}",
                               N,
                               sorted ? "sorted"sv : "unsorted"sv)};
    sumf1 = sum(flat_interp(data, flag1, flag2));
  }

  {
    test_timer_t t{std::format("flat_interp-itw fixed-cyclic {} old-N/A {}",
                               N,
                               sorted ? "sorted"sv : "unsorted"sv)};
  }

  {
    test_timer_t t{std::format("flat_interp-itw fixed-cyclic {} {}",
                               N,
                               sorted ? "sorted"sv : "unsorted"sv)};
    const auto itw = flat_interpweights(flag1, flag2);
    sumf4          = sum(flat_interp(data, flag1, flag2));
  }

  {
    test_timer_t t{std::format("flat_interp fixed-cyclic-vec {} {}",
                               N,
                               sorted ? "sorted"sv : "unsorted"sv)};
    const auto l1 =
        lagrange_interp::make_lags<N, lagrange_interp::loncross>(grid, r, 0.0);
    const auto l2 =
        lagrange_interp::make_lags<N, lagrange_interp::loncross>(grid, r2, 0.0);
    sumf2 = sum(flat_interp(data, l1, l2));
  }

  ARTS_USER_ERROR_IF(allclose(sumf0, sumf1, sumf2, sumf3, sumf4),
                     "Results of old and interpolation do not match!\n"
                     "new: {0} != {1} != {2} != {3} != {4}",
                     sumf0,
                     sumf1,
                     sumf2,
                     sumf3,
                     sumf4);
}

template <Size N>
void test_fixed_mixed(bool sorted) {
  const Vector grid{-180.0, -100.0, -50.0, 0.0, 50.0, 100.0, 179.0};

  const Matrix data =
      matpack::uniform_grid(-180, grid.size() * grid.size(), 0.1)
          .reshape(grid.size(), grid.size());

  Numeric sumf0{}, sumf1{}, sumf2{}, sumf3{}, sumf4{};

  auto r  = random_numbers(10000, -500.0, 500.0);
  auto r2 = random_numbers(r.size(), -500.0, 500.0);

  if (sorted) {
    std::sort(r.begin(), r.end());
    std::sort(r2.begin(), r2.end());
  }

  using newT1 = lagrange_interp::lag_t<N, lagrange_interp::loncross>;
  using newT2 = lagrange_interp::lag_t<N>;
  std::vector<newT1> flag1(r.size());
  std::vector<newT2> flag2(r.size());

  {
    test_timer_t t{std::format(
        "lag fixed-mixed {} {}", N, sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      flag1[i] = newT1(grid, r[i]);
      flag2[i] = newT2(grid, r2[i]);
    }
  }

  {
    test_timer_t t{std::format(
        "interp fixed-mixed {} {}", N, sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      sumf0 += interp(data, flag1[i], flag2[i]);
    }
  }

  {
    test_timer_t t{std::format(
        "interp-itw fixed-mixed {} {}", N, sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      auto itw  = interpweights(flag1[i], flag2[i]);
      sumf3    += interp(data, itw, flag1[i], flag2[i]);
    }
  }

  {
    test_timer_t t{std::format("flat_interp fixed-mixed {} {}",
                               N,
                               sorted ? "sorted"sv : "unsorted"sv)};
    sumf1 = sum(flat_interp(data, flag1, flag2));
  }

  {
    test_timer_t t{std::format("flat_interp-itw fixed-mixed {} old-N/A {}",
                               N,
                               sorted ? "sorted"sv : "unsorted"sv)};
  }

  {
    test_timer_t t{std::format("flat_interp-itw fixed-mixed {} {}",
                               N,
                               sorted ? "sorted"sv : "unsorted"sv)};
    const auto itw = flat_interpweights(flag1, flag2);
    sumf4          = sum(flat_interp(data, flag1, flag2));
  }

  {
    test_timer_t t{std::format("flat_interp fixed-mixed-vec {} {}",
                               N,
                               sorted ? "sorted"sv : "unsorted"sv)};
    const auto l1 =
        lagrange_interp::make_lags<N, lagrange_interp::loncross>(grid, r, 0.0);
    const auto l2 = lagrange_interp::make_lags<N>(grid, r2, 0.0);
    sumf2         = sum(flat_interp(data, l1, l2));
  }

  ARTS_USER_ERROR_IF(allclose(sumf0, sumf1, sumf2, sumf3, sumf4),
                     "Results of old and interpolation do not match!\n"
                     "new: {0} != {1} != {2} != {3} != {4}",
                     sumf0,
                     sumf1,
                     sumf2,
                     sumf3,
                     sumf4);
}

void test_dynamic(Index polyorder, bool sorted) {
  const Vector grid{1, 2, 3, 4, 5, 6};

  const Matrix data =
      matpack::uniform_grid(-180, grid.size() * grid.size(), 0.1)
          .reshape(grid.size(), grid.size());

  Numeric sumf0{}, sumf1{}, sumf2{}, sumf3{}, sumf4{};

  auto r  = random_numbers(10000, -5.0, 10.0);
  auto r2 = random_numbers(r.size(), -5.0, 10.0);

  if (sorted) {
    std::sort(r.begin(), r.end());
    std::sort(r2.begin(), r2.end());
  }

  using newT = lagrange_interp::lag_t<-1>;
  std::vector<newT> flag1(r.size());
  std::vector<newT> flag2(r.size());

  {
    test_timer_t t{std::format(
        "lag dynamic {} {}", polyorder, sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      flag1[i] = newT(grid, r[i], polyorder);
      flag2[i] = newT(grid, r2[i], polyorder);
    }
  }

  {
    test_timer_t t{std::format(
        "interp dynamic {} {}", polyorder, sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      sumf0 += interp(data, flag1[i], flag2[i]);
    }
  }

  {
    test_timer_t t{std::format("interp-itw dynamic {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      auto itw  = interpweights(flag1[i], flag2[i]);
      sumf3    += interp(data, itw, flag1[i], flag2[i]);
    }
  }

  {
    test_timer_t t{std::format("flat_interp dynamic {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    sumf1 = sum(flat_interp(data, flag1, flag2));
  }

  {
    test_timer_t t{std::format("flat_interp-itw dynamic {} old-N/A {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
  }

  {
    test_timer_t t{std::format("flat_interp-itw dynamic {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    const auto itw = flat_interpweights(flag1, flag2);
    sumf4          = sum(flat_interp(data, flag1, flag2));
  }

  {
    test_timer_t t{std::format("flat_interp dynamic-vec {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    const auto l1 = lagrange_interp::make_lags(grid, r, polyorder, 0.0);
    const auto l2 = lagrange_interp::make_lags(grid, r2, polyorder, 0.0);
    sumf2         = sum(flat_interp(data, l1, l2));
  }

  ARTS_USER_ERROR_IF(allclose(sumf0, sumf1, sumf2, sumf3, sumf4),
                     "Results of old and interpolation do not match!\n"
                     "new: {0} != {1} != {2} != {3} != {4}",
                     sumf0,
                     sumf1,
                     sumf2,
                     sumf3,
                     sumf4);
}

void test_dynamic_cyclic(Index polyorder, bool sorted) {
  const Vector grid{-180, -100, -50, 0, 50, 100, 179};

  const Matrix data =
      matpack::uniform_grid(-180, grid.size() * grid.size(), 0.1)
          .reshape(grid.size(), grid.size());

  Numeric sumf0{}, sumf1{}, sumf2{}, sumf3{}, sumf4{};

  auto r  = random_numbers(10000, -500.0, 500.0);
  auto r2 = random_numbers(r.size(), -500.0, 500.0);

  if (sorted) {
    std::sort(r.begin(), r.end());
    std::sort(r2.begin(), r2.end());
  }

  using newT = lagrange_interp::lag_t<-1, lagrange_interp::loncross>;
  std::vector<newT> flag1(r.size());
  std::vector<newT> flag2(r.size());

  {
    test_timer_t t{std::format("lag dynamic-cyclic {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      flag1[i] = newT(grid, r[i], polyorder);
      flag2[i] = newT(grid, r2[i], polyorder);
    }
  }

  {
    test_timer_t t{std::format("interp dynamic-cyclic {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      sumf0 += interp(data, flag1[i], flag2[i]);
    }
  }

  {
    test_timer_t t{std::format("interp-itw dynamic-cyclic {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      auto itw  = interpweights(flag1[i], flag2[i]);
      sumf3    += interp(data, itw, flag1[i], flag2[i]);
    }
  }

  {
    test_timer_t t{std::format("flat_interp dynamic-cyclic {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    sumf1 = sum(flat_interp(data, flag1, flag2));
  }

  {
    test_timer_t t{std::format("flat_interp-itw dynamic-cyclic {} old-N/A {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
  }

  {
    test_timer_t t{std::format("flat_interp-itw dynamic-cyclic {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    const auto itw = flat_interpweights(flag1, flag2);
    sumf4          = sum(flat_interp(data, flag1, flag2));
  }

  {
    test_timer_t t{std::format("flat_interp dynamic-cyclic-vec {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    const auto l1 = lagrange_interp::make_lags<lagrange_interp::loncross>(
        grid, r, polyorder, 0.0);
    const auto l2 = lagrange_interp::make_lags<lagrange_interp::loncross>(
        grid, r2, polyorder, 0.0);
    sumf2 = sum(flat_interp(data, l1, l2));
  }

  ARTS_USER_ERROR_IF(allclose(sumf0, sumf1, sumf2, sumf3, sumf4),
                     "Results of old and interpolation do not match!\n"
                     "new: {0} != {1} != {2} != {3} != {4}",
                     sumf0,
                     sumf1,
                     sumf2,
                     sumf3,
                     sumf4);
}

void test_dynamic_mixed(Index polyorder, bool sorted) {
  const Vector grid{-180, -100, -50, 0, 50, 100, 179};

  const Matrix data =
      matpack::uniform_grid(-180, grid.size() * grid.size(), 0.1)
          .reshape(grid.size(), grid.size());

  Numeric sumf0{}, sumf1{}, sumf2{}, sumf3{}, sumf4{};

  auto r  = random_numbers(10000, -500.0, 500.0);
  auto r2 = random_numbers(r.size(), -500.0, 500.0);

  if (sorted) {
    std::sort(r.begin(), r.end());
    std::sort(r2.begin(), r2.end());
  }

  using newT1 = lagrange_interp::lag_t<-1, lagrange_interp::loncross>;
  using newT2 = lagrange_interp::lag_t<-1>;
  std::vector<newT1> flag1(r.size());
  std::vector<newT2> flag2(r.size());

  {
    test_timer_t t{std::format("lag dynamic-mixed {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      flag1[i] = newT1(grid, r[i], polyorder);
      flag2[i] = newT2(grid, r2[i], polyorder);
    }
  }

  {
    test_timer_t t{std::format("interp dynamic-mixed {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      sumf0 += interp(data, flag1[i], flag2[i]);
    }
  }

  {
    test_timer_t t{std::format("interp-itw dynamic-mixed {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    for (Size i = 0; i < r.size(); ++i) {
      auto itw  = interpweights(flag1[i], flag2[i]);
      sumf3    += interp(data, itw, flag1[i], flag2[i]);
    }
  }

  {
    test_timer_t t{std::format("flat_interp dynamic-mixed {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    sumf1 = sum(flat_interp(data, flag1, flag2));
  }

  {
    test_timer_t t{std::format("flat_interp-itw dynamic-mixed {} old-N/A {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
  }

  {
    test_timer_t t{std::format("flat_interp-itw dynamic-mixed {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    const auto itw = flat_interpweights(flag1, flag2);
    sumf4          = sum(flat_interp(data, flag1, flag2));
  }

  {
    test_timer_t t{std::format("flat_interp dynamic-mixed-vec {} {}",
                               polyorder,
                               sorted ? "sorted"sv : "unsorted"sv)};
    const auto l1 = lagrange_interp::make_lags<lagrange_interp::loncross>(
        grid, r, polyorder, 0.0);
    const auto l2 = lagrange_interp::make_lags(grid, r2, polyorder, 0.0);
    sumf2         = sum(flat_interp(data, l1, l2));
  }

  ARTS_USER_ERROR_IF(allclose(sumf0, sumf1, sumf2, sumf3, sumf4),
                     "Results of old and interpolation do not match!\n"
                     "new: {0} != {1} != {2} != {3} != {4}",
                     sumf0,
                     sumf1,
                     sumf2,
                     sumf3,
                     sumf4);
}

template <Size N>
void test_reinterp_fixed() {
  const Vector grid{-180.0, -100.0, -50.0, 0.0, 50.0, 100.0, 179.0};
  const Vector ngrid1{matpack::uniform_grid(-180, 100, 1)};
  const Vector ngrid2{matpack::uniform_grid(-180, 1000, 0.4)};

  const Matrix data = matpack::uniform_grid(1, grid.size() * grid.size(), 0.1)
                          .reshape(grid.size(), grid.size());

  using NLinT = lagrange_interp::lag_t<N>;
  using NCycT = lagrange_interp::lag_t<N, lagrange_interp::loncross>;

  std::vector<NLinT> fl1;
  std::vector<NCycT> fl2;

  {
    test_timer_t t{std::format("lag reinterp fixed {} new", N)};
    {
      fl1 = lagrange_interp::make_lags<N, lagrange_interp::identity>(
          grid, ngrid1, 0.0);
      fl2 = lagrange_interp::make_lags<N, lagrange_interp::loncross>(
          grid, ngrid2, 0.0);
    }
  }

  Matrix mn, mn2;

  {
    test_timer_t t{std::format("reinterp fixed {} new", N)};
    mn = reinterp(data, fl1, fl2);
  }

  {
    test_timer_t t{std::format("reinterp-itw fixed {} new", N)};
    auto itw = reinterpweights(fl1, fl2);
    mn2      = reinterp(data, itw, fl1, fl2);
  }

  ARTS_USER_ERROR_IF(sum(mn) != sum(mn2),
                     "Reinterpolation results do not match");
}

void test_reinterp_dynamic(Index N) {
  const Vector grid{-180.0, -100.0, -50.0, 0.0, 50.0, 100.0, 179.0};
  const Vector ngrid1{matpack::uniform_grid(-180, 100, 1)};
  const Vector ngrid2{matpack::uniform_grid(-180, 1000, 0.4)};

  const Matrix data = matpack::uniform_grid(1, grid.size() * grid.size(), 0.1)
                          .reshape(grid.size(), grid.size());

  using NLinT = lagrange_interp::lag_t<-1>;
  using NCycT = lagrange_interp::lag_t<-1, lagrange_interp::loncross>;

  std::vector<NLinT> fl1;
  std::vector<NCycT> fl2;

  {
    test_timer_t t{std::format("lag reinterp dynamic {} new", N)};
    {
      fl1 = lagrange_interp::make_lags<lagrange_interp::identity>(
          grid, ngrid1, N, 0.0);
      fl2 = lagrange_interp::make_lags<lagrange_interp::loncross>(
          grid, ngrid2, N, 0.0);
    }
  }

  Matrix mn, mn2;

  {
    test_timer_t t{std::format("reinterp dynamic {} new", N)};
    mn = reinterp(data, fl1, fl2);
  }

  {
    test_timer_t t{std::format("reinterp-itw dynamic {} new", N)};
    auto itw = reinterpweights(fl1, fl2);
    mn2      = reinterp(data, itw, fl1, fl2);
  }

  ARTS_USER_ERROR_IF(sum(mn) != sum(mn2),
                     "Reinterpolation results do not match");
}
}  // namespace

int main() {
  for (Index i = 0; i < 5; ++i) {
    test_fixed<0>(true);
    test_fixed<1>(true);
    test_fixed<2>(true);
    test_fixed<3>(true);
    test_fixed<4>(true);
    test_fixed<5>(true);

    test_fixed<0>(false);
    test_fixed<1>(false);
    test_fixed<2>(false);
    test_fixed<3>(false);
    test_fixed<4>(false);
    test_fixed<5>(false);

    test_dynamic(0, true);
    test_dynamic(1, true);
    test_dynamic(2, true);
    test_dynamic(3, true);
    test_dynamic(4, true);
    test_dynamic(5, true);

    test_dynamic(0, false);
    test_dynamic(1, false);
    test_dynamic(2, false);
    test_dynamic(3, false);
    test_dynamic(4, false);
    test_dynamic(5, false);

    test_fixed_cyclic<0>(true);
    test_fixed_cyclic<1>(true);
    test_fixed_cyclic<2>(true);
    test_fixed_cyclic<3>(true);
    test_fixed_cyclic<4>(true);
    test_fixed_cyclic<5>(true);
    test_fixed_cyclic<6>(true);

    test_fixed_cyclic<0>(false);
    test_fixed_cyclic<1>(false);
    test_fixed_cyclic<2>(false);
    test_fixed_cyclic<3>(false);
    test_fixed_cyclic<4>(false);
    test_fixed_cyclic<5>(false);
    test_fixed_cyclic<6>(false);

    test_fixed_mixed<0>(true);
    test_fixed_mixed<1>(true);
    test_fixed_mixed<2>(true);
    test_fixed_mixed<3>(true);
    test_fixed_mixed<4>(true);
    test_fixed_mixed<5>(true);
    test_fixed_mixed<6>(true);

    test_fixed_mixed<0>(false);
    test_fixed_mixed<1>(false);
    test_fixed_mixed<2>(false);
    test_fixed_mixed<3>(false);
    test_fixed_mixed<4>(false);
    test_fixed_mixed<5>(false);
    test_fixed_mixed<6>(false);

    test_dynamic_cyclic(0, false);
    test_dynamic_cyclic(1, false);
    test_dynamic_cyclic(2, false);
    test_dynamic_cyclic(3, false);
    test_dynamic_cyclic(4, false);
    test_dynamic_cyclic(5, false);
    test_dynamic_cyclic(6, false);

    test_dynamic_cyclic(0, true);
    test_dynamic_cyclic(1, true);
    test_dynamic_cyclic(2, true);
    test_dynamic_cyclic(3, true);
    test_dynamic_cyclic(4, true);
    test_dynamic_cyclic(5, true);
    test_dynamic_cyclic(6, true);

    test_dynamic_mixed(0, false);
    test_dynamic_mixed(1, false);
    test_dynamic_mixed(2, false);
    test_dynamic_mixed(3, false);
    test_dynamic_mixed(4, false);
    test_dynamic_mixed(5, false);
    test_dynamic_mixed(6, false);

    test_dynamic_mixed(0, true);
    test_dynamic_mixed(1, true);
    test_dynamic_mixed(2, true);
    test_dynamic_mixed(3, true);
    test_dynamic_mixed(4, true);
    test_dynamic_mixed(5, true);
    test_dynamic_mixed(6, true);

    test_reinterp_fixed<0>();
    test_reinterp_fixed<1>();
    test_reinterp_fixed<2>();
    test_reinterp_fixed<3>();
    test_reinterp_fixed<4>();
    test_reinterp_fixed<5>();

    test_reinterp_dynamic(0);
    test_reinterp_dynamic(1);
    test_reinterp_dynamic(2);
    test_reinterp_dynamic(3);
    test_reinterp_dynamic(4);
    test_reinterp_dynamic(5);
  }

  print_time_points();
}
