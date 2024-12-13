#pragma once

#include <disort.h>
#include <matpack.h>

#include <iostream>

#include "matpack_iter.h"
#include "nonstd.h"

inline bool is_good(const auto& a, const auto& b) {
  if (a.shape() != b.shape()) {
    std::cerr << "!!!\n\nBad shapes\n\n!!!\n";
    return false;
  }

  return std::transform_reduce(
      a.elem_begin(),
      a.elem_end(),
      b.elem_begin(),
      true,
      [](bool first, bool second) { return first and second; },
      [](Numeric first, Numeric second) {
        if (nonstd::isnan(first) or nonstd::isnan(second)) return false;
        if (first == 0.0 and second == 0.0) return true;
        const Numeric ratio = std::abs(first / second - 1);
        return ratio < 1e-6;
      });
}

inline Tensor3 compute_u(const disort::main_data& dis,
                         const Vector& taus,
                         const Vector& phis,
                         const bool nt_corr) {
  Tensor3 u(phis.size(), taus.size(), dis.mu().size());
  disort::u_data u_data;
  Vector ims;
  disort::tms_data tms_data;

  for (Index j = 0; j < taus.size(); j++) {
    for (Index i = 0; i < phis.size(); i++) {
      if (nt_corr) {
        dis.u_corr(u_data, ims, tms_data, taus[j], phis[i]);
      } else {
        dis.u(u_data, taus[j], phis[i]);
      }
      u[i, j, joker] = u_data.intensities;
    }
  }
  return u;
}

inline Matrix compute_u0(const disort::main_data& dis, const Vector& taus) {
  Matrix u0(taus.size(), dis.mu().size());
  disort::u0_data u0_data;

  for (Index j = 0; j < taus.size(); j++) {
    dis.u0(u0_data, taus[j]);
    u0[j] = u0_data.u0;
  }
  return u0;
}

inline std::tuple<Vector, Vector, Vector> compute_flux(
    const disort::main_data& dis, const Vector& taus) {
  Vector flux_up(taus.size()), flux_down_diffuse(taus.size()),
      flux_down_direct(taus.size());
  disort::flux_data flux_data;

  for (Index j = 0; j < taus.size(); j++) {
    auto [ds, dr]        = dis.flux_down(flux_data, taus[j]);
    flux_up[j]           = dis.flux_up(flux_data, taus[j]);
    flux_down_diffuse[j] = ds;
    flux_down_direct[j]  = dr;
  }
  return {flux_up, flux_down_diffuse, flux_down_direct};
}

inline void compare(const std::string_view name,
                    const disort::main_data& dis,
                    const Vector& taus,
                    const Vector& phis,
                    const Tensor3& u,
                    const Matrix& u0,
                    const Vector& flux_down_diffuse,
                    const Vector& flux_down_direct,
                    const Vector& flux_up,
                    const bool nt_corr) {
  const auto u_arts  = compute_u(dis, taus, phis, nt_corr);
  const auto u0_arts = compute_u0(dis, taus);
  const auto [flux_up_arts, flux_down_diffuse_arts, flux_down_direct_arts] =
      compute_flux(dis, taus);

  std::cout << std::format(
      "u_arts:\n{:B,}\nu0_arts:\n{:B,}\nflux_up_arts:\n{:B,}\nflux_down_diffuse_arts:\n{:B,}\nflux_down_direct_arts:\n{:B,}\n",
      u_arts.flat_view(),
      u0_arts.flat_view(),
      flux_up_arts,
      flux_down_diffuse_arts,
      flux_down_direct_arts);

  ARTS_USER_ERROR_IF(not is_good(u_arts, u), "Failed u in test {}", name);
  ARTS_USER_ERROR_IF(not is_good(u0_arts, u0), "Failed u0 in test {}", name);
  ARTS_USER_ERROR_IF(
      not is_good(flux_up_arts, flux_up), "Failed flux_up in test {}", name);
  ARTS_USER_ERROR_IF(not is_good(flux_down_diffuse_arts, flux_down_diffuse),
                     "Failed flux_down_diffuse in test {}",
                     name);
  ARTS_USER_ERROR_IF(not is_good(flux_down_direct_arts, flux_down_direct),
                     "Failed flux_down_direct in test {}",
                     name);
};

inline void flat_print(const auto& a, const auto& b) {
  ARTS_USER_ERROR_IF(a.shape() != b.shape(), "Failed shape comparison");
  const auto av = a.flat_view();
  const auto bv = b.flat_view();

  for (Index i = 0; i < a.size(); i++) {
    std::cout << i << ' ' << av[i] << ' ' << bv[i] << ' ' << (bv[i] / av[i] - 1)
              << '\n';
  }
}