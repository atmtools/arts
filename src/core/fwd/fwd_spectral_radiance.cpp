#include "fwd_spectral_radiance.h"

#include <path_point.h>
#include <physics_funcs.h>

#include <algorithm>
#include <memory>
#include <ostream>

#include "arts_constants.h"
#include "arts_omp.h"
#include "debug.h"
#include "fwd_propmat.h"
#include "rtepack.h"

namespace fwd {
spectral_radiance_1d::spectral_radiance_1d(
    AscendingGrid alt_,
    Numeric lat,
    Numeric lon,
    const AtmField& atm_,
    const SurfaceField& surf,
    std::shared_ptr<AbsorptionBands> lines,
    std::shared_ptr<ArrayOfCIARecord> cia,
    std::shared_ptr<ArrayOfXsecRecord> xsec,
    std::shared_ptr<PredefinedModelData> predef,
    Numeric ciaextrap,
    Index ciarobust)
    : alt(std::move(alt_)),
      latitude(lat),
      longitude(lon),
      atm(alt.size()),
      spectral_radiance_surface(
          [surf = surf.at(lat, lon)](Numeric f, Vector2) -> Stokvec {
            return planck(f, surf.temperature);
          }),
      spectral_radiance_space([](Numeric f) -> Stokvec {
        return planck(f, Constant::cosmic_microwave_background_temperature);
      }) {
  ARTS_USER_ERROR_IF(alt.size() == 0, "Must have a sized atmosphere")

  if (arts_omp_in_parallel()) {
    for (Index i = 0; i < alt.size(); i++) {
      atm[i] = std::make_shared<AtmPoint>(atm_.at(alt[i], lat, lon));
    }
  } else {
    String error{};
#pragma omp parallel for
    for (Index i = 0; i < alt.size(); i++) {
      try {
        atm[i] = std::make_shared<AtmPoint>(atm_.at(alt[i], lat, lon));
      } catch (std::exception& e) {
#pragma omp critical
        error += e.what() + String("\n");
      }
    }
    ARTS_USER_ERROR_IF(not error.empty(), error)
  }

  pm.resize(alt.size(),
            propmat(atm[0],
                    std::move(lines),
                    std::move(cia),
                    std::move(xsec),
                    std::move(predef),
                    ciaextrap,
                    ciarobust));

  if (arts_omp_in_parallel()) {
    for (Index i = 1; i < alt.size(); i++) {
      pm[i].set_atm(atm[i]);
    }
  } else {
    String error{};
#pragma omp parallel for
    for (Index i = 1; i < alt.size(); i++) {
      try {
        pm[i].set_atm(atm[i]);
      } catch (std::exception& e) {
#pragma omp critical
        error += e.what() + String("\n");
      }
    }
    ARTS_USER_ERROR_IF(not error.empty(), error)
  }
}

Stokvec spectral_radiance_1d::operator()(const Numeric frequency,
                                         const Vector2 los,
                                         const PathVectorView path) const {
  ARTS_ASSERT(path.size() > 0)

  const auto B = [f = frequency, this](const Size i) -> Stokvec {
    return {planck(f, atm[i]->temperature), 0, 0, 0};
  };

  const auto Iback = [f = frequency, l = los, this](bool surface) -> Stokvec {
    return surface ? spectral_radiance_surface(f, l)
                   : spectral_radiance_space(f);
  };

  if (path.front().i >= atm.size()) {
    return Iback(path.front().r < 0.0);
  }

  auto [K, N] = pm[path.front().i](frequency, los);
  const Stokvec J = inv(K) * N + B(path.front().i);
  Muelmat T = path.front().r == 0.0 ? 1.0 : exp(K, path.front().r);
  Stokvec I = (1.0 - T) * J;

  for (auto& [i, r] : path | std::views::drop(1)) {
    if (i >= atm.size()) {
      I += T * Iback(path.front().r < 0.0);
      return I;
    }

    const auto [Ki, Ni] = pm[i](frequency, los);
    const Stokvec Ji = inv(Ki) * Ni + B(i);
    const Muelmat Ti = T * exp(avg(Ki, K), r);

    if (Ti(0, 0) < 1e-6) {
      I += Ti * Ji;
      return I;
    }

    I += (T - Ti) * Ji;

    K = Ki;
    T = Ti;
  }

  return I;
}

void spectral_radiance_1d::operator()(StokvecVectorView out,
                                      const Numeric frequency,
                                      const Vector2 los,
                                      const PathVectorView path) const {
  ARTS_ASSERT(path.size() > 0)
  ARTS_ASSERT(out.size() == alt.size())

  const auto B = [f = frequency, this](const Size i) -> Stokvec {
    return {planck(f, atm[i]->temperature), 0, 0, 0};
  };

  const auto Iback = [f = frequency, l = los, this](bool surface) -> Stokvec {
    return surface ? spectral_radiance_surface(f, l)
                   : spectral_radiance_space(f);
  };

  if (path.front().i >= atm.size()) {
    out.back() = Iback(path.front().r < 0.0);
  }

  auto* ptr = std::ranges::find_if(
      path,
      [n = atm.size()](const Size& is) { return is >= n; },
      &PosDistance::i);

  const bool up_looking = ptr->r > 0.0;
  Stokvec I = Iback(up_looking);

  --ptr;

  auto [K, N] = pm[ptr->i](frequency, los);
  Stokvec J = inv(K) * N + B(ptr->i);
  Muelmat T = ptr->r == 0.0 ? 1.0 : exp(K, ptr->r);

  I = T * (I - J) + J;
  out[ptr->i] = I;

  while (--ptr >= path.begin()) {
    const auto [i, r] = *ptr;

    const auto [Ki, Ni] = pm[i](frequency, los);
    const auto Ji = inv(Ki) * Ni + B(i);
    const auto Ti = exp(avg(Ki, K), r);

    I = Ti * (I - avg(J, Ji)) + avg(J, Ji);
    out[i] = I;

    K = Ki;
    J = Ji;
  }
}

std::ostream& operator<<(std::ostream& os,
                         const spectral_radiance_1d::PosDistance& sr) {
  return os << "i: " << sr.i << ", r: " << sr.r << "\n";
}

std::ostream& operator<<(std::ostream& os, const spectral_radiance_1d& sr) {
  return os << "Spectral radiance operator:\n"
            << "  Altitude grid: " << sr.alt << "\n";
}

spectral_radiance_1d::PathVector spectral_radiance_1d::geometric_planar(
    const Numeric altitude, const Numeric zenith) const {
  ARTS_USER_ERROR_IF(zenith == 90.0,
                     "Not a valid zenith angle for geometric planar radiation")
  ARTS_USER_ERROR_IF(alt.size() == 0, "No altitude grid")
  ARTS_USER_ERROR_IF(alt[0] > altitude, "Subsurface altitude")

  std::vector<PosDistance> path;
  path.reserve(alt.size());

  const bool up_looking = zenith > 90.0;
  const Numeric csc = std::abs(1.0 / Conversion::cosd(zenith));

  Size i0 = std::distance(
      alt.begin(), std::find_if(alt.begin(), alt.end(), [altitude](Numeric a) {
        return a >= altitude;
      }));

  if (i0 >= static_cast<Size>(alt.size())) {
    path.emplace_back(alt.size(), up_looking ? 1.0 : -1.0);
    return path;
  }

  i0 -= not up_looking and altitude < alt[i0];

  path.emplace_back(i0, std::abs(alt[i0] - altitude) * csc);

  if (up_looking) {
    while (path.back().i < static_cast<Size>(alt.size() - 1)) {
      auto x = path.back();
      path.emplace_back(x.i + 1, std::abs(alt[x.i + 1] - alt[x.i]) * csc);
    }
    path.emplace_back(alt.size(), 1.0);
  } else {
    while (path.back().i > 0) {
      auto x = path.back();
      path.emplace_back(x.i - 1, std::abs(alt[x.i - 1] - alt[x.i]) * csc);
    }
    path.emplace_back(alt.size(), -1.0);
  }

  return std::move(path);
}

ArrayOfAtmPoint spectral_radiance_1d::get_atm(const PathVectorView path) const {
  ArrayOfAtmPoint out;
  for (const auto& pd : path) {
    if (pd.i >= atm.size()) break;
    out.push_back(*atm[pd.i]);
  }
  return out;
}

ExhaustiveConstVectorView spectral_radiance_1d::altitude_grid() const {
  return alt;
}

ExhaustiveConstVectorView spectral_radiance_1d::latitude_grid() const {
  return ExhaustiveConstVectorView{latitude};
}

ExhaustiveConstVectorView spectral_radiance_1d::longitude_grid() const {
  return ExhaustiveConstVectorView{longitude};
}
}  // namespace fwd