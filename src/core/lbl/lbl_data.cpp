#include "lbl_data.h"

#include <cmath>

#include "arts_constants.h"
#include "arts_constexpr_math.h"
#include "quantum_numbers.h"

//! In CPP file
using Constant::c;
using Constant::k;
using Constant::pi;
using Math::pow2;
using Math::pow3;
using Math::pow4;
using std::exp;

namespace lbl {
Numeric line::s(Numeric T, Numeric Q) const noexcept {
  return a * pow2(c) * gl * exp(-e0 / (k * T)) / (8 * pi * pow3(f0) * Q);
}

Numeric line::ds_dT(Numeric T, Numeric Q, Numeric dQ_dT) const noexcept {
  return a * pow2(c) * gl * (e0 * Q - k * pow2(T) * dQ_dT) *
         exp(-e0 / (k * T)) / (8 * pi * pow3(f0) * k * pow2(T) * pow2(Q));
}

Numeric line::ds_de0(Numeric T, Numeric Q) const noexcept {
  return -a * pow2(c) * gl * exp(-e0 / (k * T)) /
         (8 * pi * pow3(f0) * k * T * Q);
}

Numeric line::ds_df0(Numeric T, Numeric Q) const noexcept {
  return -3 * a * pow2(c) * gl * exp(-e0 / (k * T)) / (8 * pi * pow4(f0) * Q);
}

Numeric line::ds_da(Numeric T, Numeric Q) const noexcept {
  return pow2(c) * gl * exp(-e0 / (k * T)) / (8 * pi * pow3(f0) * Q);
}
}  // namespace lbl
