#include "fwd_spectral_radiance.h"

#include <path_point.h>
#include <physics_funcs.h>

#include <iostream>
#include <memory>

#include "arts_constants.h"
#include "arts_omp.h"
#include "debug.h"
#include "fwd_propmat.h"
#include "rtepack.h"

namespace fwd {
spectral_radiance::spectral_radiance(
    AscendingGrid alt_,
    Numeric lat,
    Numeric lon,
    const AtmField& atm_,
    const SurfaceField& surf_,
    std::shared_ptr<AbsorptionBands> lines,
    std::shared_ptr<ArrayOfCIARecord> cia,
    std::shared_ptr<ArrayOfXsecRecord> xsec,
    std::shared_ptr<PredefinedModelData> predef,
    Numeric ciaextrap,
    Index ciarobust)
    : alt(std::move(alt_)), atm(alt.size()), surf{surf_.at(lat, lon)} {
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

spectral_radiance::PosDistance closeby_altitude_position(const Vector& alt,
                                                         const Numeric z) {
  const auto p = std::min(std::lower_bound(alt.begin(), alt.end(), z), alt.end() - 1);
  return {.i = static_cast<Size>(std::distance(alt.begin(), p)),
          .r = *p - z};
}

Stokvec spectral_radiance::operator()(const Numeric frequency,
                                      const Vector2 los,
                                      const PathVectorView path) const {
  ARTS_ASSERT(path.size() > 0)

  constexpr Numeric Tsp = Constant::cosmic_microwave_background_temperature;

  const auto B = [f = frequency, this](const Size i) {
    return Stokvec{planck(f, atm[i]->temperature), 0, 0, 0};
  };
  const auto Iback = [f = frequency](Numeric t) {
    return Stokvec{planck(f, t), 0, 0, 0};
  };

  auto [K, N] = pm[path.front().i](frequency, los);
  const Stokvec J = inv(K) * N + B(path.front().i);
  Muelmat T = path.front().r == 0.0 ? 1.0 : exp(K, path.front().r);
  Stokvec I = (1.0 - T) * J;

  for (auto& [i, r] : path | std::views::drop(1)) {
    if (i >= atm.size()) {
      I += T * Iback(r < 0.0 ? surf.temperature : Tsp);
      return I;
    }

    const auto [Ki, Ni] = pm[i](frequency, los);
    const Stokvec Ji = inv(Ki) * Ni + B(i);
    const Muelmat Ti = exp(avg(Ki, K), r) * T;

    if (Ti(0, 0) <= 1e-6) {
      I += T * Ji;
      return I;
    }

    I += (T - Ti) * Ji;

    K = Ki;
    T = Ti;
  }

  return I;
}

std::ostream& operator<<(std::ostream& os,
                         const spectral_radiance::PosDistance& sr) {
  return os << "i: " << sr.i << ", r: " << sr.r << "\n";
}

std::ostream& operator<<(std::ostream& os, const spectral_radiance& sr) {
  return os << "Spectral radiance operator:\n"
            << "  Altitude grid: " << sr.alt << "\n";
}

spectral_radiance::PathVector spectral_radiance::geometric_planar(
    const Numeric altitude, const Numeric zenith) const {
  ARTS_USER_ERROR_IF(zenith == 90.0,
                     "Not a valid zenith angle for geometric planar radiation")

  std::vector<PosDistance> path;
  path.reserve(alt.size());

  const bool up_looking = zenith > 90.0;
  const Numeric csc = std::abs(1.0 / std::cos(zenith));

  if (const auto start_posd = closeby_altitude_position(alt, altitude);
      up_looking) {
    path.emplace_back(start_posd.i + 1,
                      std::abs(altitude - alt[start_posd.i + 1]) * csc);
  } else {
    path.emplace_back(start_posd.i, std::abs(start_posd.r) * csc);
  }

  if (up_looking) {
    while (path.back().i < static_cast<Size>(alt.size()) - 1) {
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
}  // namespace fwd