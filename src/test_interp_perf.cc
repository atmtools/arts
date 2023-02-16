#include "artstime.h"
#include "interp.h"
#include "interpolation.h"
#include "matpack_math.h"
#include "matpack_math_extra.h"

#include <concepts>
#include <cstdlib>
#include <ostream>
#include <stdexcept>

using LinearRuntime = my_interp::Lagrange<-1>;
using LinearFixed1 = my_interp::Lagrange<1>;
using LinearFixed5 = my_interp::Lagrange<5>;

auto linear_lag_runtime_raw(const Vector& xs, const Vector& xi, Index order) {
  return my_interp::lagrange_interpolation_list<LinearRuntime>(xs, xi, order);
}

auto linear_lag_runtime_view(const ConstVectorView& xs, const ConstVectorView& xi, Index order) {
  return my_interp::lagrange_interpolation_list<LinearRuntime>(xs, xi, order);
}

auto linear_lag_compile_raw_1lin(const Vector& xs, const Vector& xi) {
  return my_interp::lagrange_interpolation_list<LinearFixed1>(xs, xi);
}

auto linear_lag_compile_view_1lin(const ConstVectorView& xs, const ConstVectorView& xi) {
  return my_interp::lagrange_interpolation_list<LinearFixed1>(xs, xi);
}

auto linear_lag_compile_raw_5lin(const Vector& xs, const Vector& xi) {
  return my_interp::lagrange_interpolation_list<LinearFixed5>(xs, xi);
}

auto linear_lag_compile_view_5lin(const ConstVectorView& xs, const ConstVectorView& xi) {
  return my_interp::lagrange_interpolation_list<LinearFixed5>(xs, xi);
}

template <typename ... SomeLags>
auto reinterp_direct(const auto& f, const SomeLags&...lag) {
  return reinterp(f, lag...);
}

template <typename ... SomeLags>
auto reinterp_indirect(const auto& f, const SomeLags&...lag) {
  auto iw = interpweights(lag...);
  return reinterp(f, iw, lag...);
}

struct Timing {
  std::string_view name;
  Timing(const char * c) : name(c) {}
  TimeStep dt{};
  template <typename Function> void operator()(Function&& f) {
    Time start{};
    f();
    Time end{};
    dt = end - start;
  }
};

std::ostream& operator<<(std::ostream& os, const std::vector<Timing>& vt) {
  for (auto& t: vt) if (t.name not_eq "dummy") os << t.name  << " : " << t.dt << '\n';
  return os;
}

std::vector<Timing> test_2x_reinterp_vector(Index n) {
  Vector xi=uniform_grid(0, n, 1.0/static_cast<Numeric>(n-1));
  Vector yi(n); transform(yi, sin, xi);
  Vector xs=uniform_grid(0, 2*n-1, 1.0/static_cast<Numeric>(2*n-2));
  
  std::vector<Numeric> xvec;
  std::vector<Timing> out;
  Numeric X;

  //! Startup
  out.emplace_back("dummy")([&](){
    X = sum(xi) + sum(xs) + xi*xi + xs*xs;
  });
  xvec.push_back(X);
  
  //// 1D LINEAR

  out.emplace_back("runtime-order-1-raw-indirect")([&](){
    auto lag = linear_lag_runtime_raw(xs, xi, 1);
    auto ys = reinterp_indirect(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("runtime-order-1-view-indirect")([&](){
    auto lag = linear_lag_runtime_view(xs, xi, 1);
    auto ys = reinterp_indirect(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("compile-order-1-raw-indirect")([&](){
    auto lag = linear_lag_compile_raw_1lin(xs, xi);
    auto ys = reinterp_indirect(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("compile-order-1-view-indirect")([&](){
    auto lag = linear_lag_compile_view_1lin(xs, xi);
    auto ys = reinterp_indirect(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("runtime-order-1-raw-direct")([&](){
    auto lag = linear_lag_runtime_raw(xs, xi, 1);
    auto ys = reinterp_direct(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("runtime-order-1-view-direct")([&](){
    auto lag = linear_lag_runtime_view(xs, xi, 1);
    auto ys = reinterp_direct(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("compile-order-1-raw-direct")([&](){
    auto lag = linear_lag_compile_raw_1lin(xs, xi);
    auto ys = reinterp_direct(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("compile-order-1-view-direct")([&](){
    auto lag = linear_lag_compile_view_1lin(xs, xi);
    auto ys = reinterp_direct(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("gridpos")([&](){
    ArrayOfGridPos gp(xs.size());
    gridpos(gp, xi, xs);
    Matrix iw(gp.size(), 2);
    interpweights(iw, gp);
    Vector ys(xs.size());
    interp(ys, iw, yi, gp);
    X = sum(ys);
  });
  xvec.push_back(X);

  // 1D 5order

  out.emplace_back("runtime-order-5-raw-indirect")([&](){
    auto lag = linear_lag_runtime_raw(xs, xi, 5);
    auto ys = reinterp_indirect(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("runtime-order-5-view-indirect")([&](){
    auto lag = linear_lag_runtime_view(xs, xi, 5);
    auto ys = reinterp_indirect(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("compile-order-5-raw-indirect")([&](){
    auto lag = linear_lag_compile_raw_5lin(xs, xi);
    auto ys = reinterp_indirect(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("compile-order-5-view-indirect")([&](){
    auto lag = linear_lag_compile_view_5lin(xs, xi);
    auto ys = reinterp_indirect(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("runtime-order-5-raw-direct")([&](){
    auto lag = linear_lag_runtime_raw(xs, xi, 5);
    auto ys = reinterp_direct(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("runtime-order-5-view-direct")([&](){
    auto lag = linear_lag_runtime_view(xs, xi, 5);
    auto ys = reinterp_direct(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("compile-order-5-raw-direct")([&](){
    auto lag = linear_lag_compile_raw_5lin(xs, xi);
    auto ys = reinterp_direct(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  out.emplace_back("compile-order-5-view-direct")([&](){
    auto lag = linear_lag_compile_view_5lin(xs, xi);
    auto ys = reinterp_direct(yi, lag);
    X = sum(ys);
  });
  xvec.push_back(X);

  return out;
}

std::vector<Timing> test_linear_startup_cost(Index n) {
  Vector xi=uniform_grid(0, n, 1.0/static_cast<Numeric>(n-1));
  Vector yi(n); transform(yi, sin, xi);
  Vector xs=uniform_grid(0, 2*n-1, 1.0/static_cast<Numeric>(2*n-2));
  
  std::vector<Numeric> xvec;
  std::vector<Timing> out;
  Numeric X;
    //! Startup
  out.emplace_back("dummy")([&](){
    X = sum(xi) + sum(xs) + xi*xi + xs*xs;
  });
  xvec.push_back(X);

  out.emplace_back("gridpos")([&](){
    ArrayOfGridPos gp(xs.size());
    gridpos(gp, xi, xs);
    X = std::transform_reduce(gp.begin(), gp.end(), Numeric{}, std::plus<>(), [](auto& g){return g.fd[0];});
  });
  xvec.push_back(X);

  out.emplace_back("compile-order-1")([&](){
    const auto lag = my_interp::lagrange_interpolation_list<LinearFixed1>(xs, xi);
    X = std::transform_reduce(lag.begin(), lag.end(), Numeric{}, std::plus<>(), [](auto& g){return g.lx[0];});
  });
  xvec.push_back(X);

  out.emplace_back("runtime-order-1")([&](){
    const auto lag = linear_lag_runtime_raw(xs, xi, 1);
    X = std::transform_reduce(lag.begin(), lag.end(), Numeric{}, std::plus<>(), [](auto& g){return g.lx[0];});
  });
  xvec.push_back(X);

  return out;
}

std::vector<Timing> test_linear_interpweights_cost(Index n) {
  Vector xi=uniform_grid(0, n, 1.0/static_cast<Numeric>(n-1));
  Vector yi(n); transform(yi, sin, xi);
  Vector xs=uniform_grid(0, 2*n-1, 1.0/static_cast<Numeric>(2*n-2));
  
  std::vector<Numeric> xvec;
  std::vector<Timing> out;
  Numeric X;
    //! Startup
  out.emplace_back("dummy")([&](){
    X = sum(xi) + sum(xs) + xi*xi + xs*xs;
  });
  xvec.push_back(X);

  ArrayOfGridPos gp(xs.size());
  gridpos(gp, xi, xs);
  const auto lagc = my_interp::lagrange_interpolation_list<LinearFixed1>(xs, xi);
  const auto lagr = linear_lag_runtime_raw(xs, xi, 1);

  out.emplace_back("gridpos")([&](){
    Matrix iw(gp.size(), 2);
    interpweights(iw, gp);
    X = iw(0, 0);
  });
  xvec.push_back(X);

  out.emplace_back("compile-order-1")([&]() {
    const auto iw = interpweights(lagc);
    X = iw[0][0];
  });
  xvec.push_back(X);

  out.emplace_back("runtime-order-1")([&](){
    const auto iw = interpweights(lagr);
    X = iw[0][0];
  });
  xvec.push_back(X);

  return out;
}

int main(int argc, char** c) {
    std::array <Index, 3> N;
  if (static_cast<std::size_t>(argc) < 1 + 1 + N.size()) {
    std::cerr << "Expects PROGNAME NREPEAT NSIZE..., wehere NSIZE is " << N.size() << " indices\n";
    return EXIT_FAILURE;
  }

  const auto n = static_cast<Index>(std::atoll(c[1]));
  for (std::size_t i=0; i<N.size(); i++)  N[i] = static_cast<Index>(std::atoll(c[2 + i]));

  for (Index i=0; i<n; i++) {
    std::cout << N[0] << " input test_2x_reinterp_vector\n" << test_2x_reinterp_vector(N[0]) << '\n';
    std::cout << N[1] << " input test_linear_startup_cost\n" << test_linear_startup_cost(N[1]) << '\n';
    std::cout << N[2] << " input test_linear_interpweights_cost\n" << test_linear_interpweights_cost(N[2]) << '\n';
  }
}
