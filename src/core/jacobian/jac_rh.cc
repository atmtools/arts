#include "jac_rh.h"

#include <debug.h>

#include <stdexcept>

#include "jacobian.h"
#include "operators.h"

Vector rhfwd::operator()(ConstVectorView x, const AtmField& atm) const try {
  const auto& t = atm[AtmKey::t];
  const auto& p = atm[AtmKey::p];

  const auto& [alt, lat, lon] = atm[self_key].get<GeodeticField3>().grids;

  ARTS_USER_ERROR_IF((alt.size() * lat.size() * lon.size()) != x.size(),
                     R"(Mismatch expected size:

x must have the size of alt x lat x lon

x.size():   {}
alt.size(): {}
lat.size(): {}
lon.size(): {}
)",
                     x.size(),
                     alt.size(),
                     lat.size(),
                     lon.size())

  Vector out{x};

  Size i = 0;
  for (auto& a : alt) {
    for (auto& la : lat) {
      for (auto& lo : lon) {
        out[i] *= p.at(a, la, lo) / psat(t.at(a, la, lo));
        if (self_fix) out[i] = std::max(out[i], 0.0);
        i++;
      }
    }
  }

  return out;
} catch (std::invalid_argument& e) {
  throw std::runtime_error(
      std::format("Error in rhfwd for key {}:\n{}", self_key, e.what()));
} catch (...) {
  throw;
}

Vector rhinv::operator()(ConstVectorView x, const AtmField& atm) const try {
  const auto& t = atm[AtmKey::t];
  const auto& p = atm[AtmKey::p];

  const auto& [alt, lat, lon] = atm[self_key].get<GeodeticField3>().grids;

  ARTS_USER_ERROR_IF((alt.size() * lat.size() * lon.size()) != x.size(),
                     R"(Mismatch expected size:

x must have the size of alt x lat x lon

x.size():   {}
alt.size(): {}
lat.size(): {}
lon.size(): {}
)",
                     x.size(),
                     alt.size(),
                     lat.size(),
                     lon.size())

  Vector out{x};

  Size i = 0;
  for (auto& a : alt) {
    for (auto& la : lat) {
      for (auto& lo : lon) {
        out[i] *= psat(t.at(a, la, lo)) / p.at(a, la, lo);
        if (self_fix) out[i] = std::max(out[i], 0.0);
        i++;
      }
    }
  }

  return out;
} catch (std::invalid_argument& e) {
  throw std::runtime_error(
      std::format("Error in rhfwd for key {}:\n{}", self_key, e.what()));
} catch (...) {
  throw;
}

Matrix rhinv::operator()(ConstMatrixView dy,
                         ConstVectorView x,
                         const AtmField& atm) const try {
  const auto& t = atm[AtmKey::t];
  const auto& p = atm[AtmKey::p];

  const auto& [alt, lat, lon] = atm[self_key].get<GeodeticField3>().grids;

  ARTS_USER_ERROR_IF((alt.size() * lat.size() * lon.size()) != x.size() or
                         static_cast<Size>(dy.ncols()) != x.size(),
                     R"(Mismatch expected size:

x must have the size of alt x lat x lon, the same as the column size of dy

x.size():   {}
alt.size(): {}
lat.size(): {}
lon.size(): {}
dy.shape(): {:B,}
)",
                     x.size(),
                     alt.size(),
                     lat.size(),
                     lon.size(),
                     dy.shape())

  Matrix out{dy};
  Size i = 0;

  for (auto& a : alt) {
    for (auto& la : lat) {
      for (auto& lo : lon) {
        out[joker, i] *= psat(t.at(a, la, lo)) / p.at(a, la, lo);
        i++;
      }
    }
  }

  return out;
} catch (std::invalid_argument& e) {
  throw std::runtime_error(
      std::format("Error in rhfwd for key {}:\n{}", self_key, e.what()));
} catch (...) {
  throw;
}

void make_rhfit(Jacobian::AtmTarget& x,
                const AtmField&,
                const NumericUnaryOperator& unary_op,
                const bool self_fix) {
  const rhfwd rfwd{.psat = unary_op, .self_key = x.type, .self_fix = self_fix};
  const rhinv rinv{.psat = unary_op, .self_key = x.type, .self_fix = self_fix};
  x.inverse_state    = rinv;
  x.inverse_jacobian = rinv;
  x.transform_state  = rfwd;
}
