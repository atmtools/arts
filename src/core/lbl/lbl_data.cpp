#include "lbl_data.h"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <limits>

#include "arts_constants.h"
#include "arts_constexpr_math.h"
#include "double_imanip.h"
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
  return a * pow2(c) * gu * exp(-e0 / (k * T)) / (8 * pi * pow3(f0) * Q);
}

Numeric line::ds_dT(Numeric T, Numeric Q, Numeric dQ_dT) const noexcept {
  return a * pow2(c) * gu * (e0 * Q - k * pow2(T) * dQ_dT) *
         exp(-e0 / (k * T)) / (8 * pi * pow3(f0) * k * pow2(T) * pow2(Q));
}

Numeric line::ds_de0(Numeric T, Numeric Q) const noexcept {
  return -a * pow2(c) * gu * exp(-e0 / (k * T)) /
         (8 * pi * pow3(f0) * k * T * Q);
}

Numeric line::ds_df0(Numeric T, Numeric Q) const noexcept {
  return -3 * a * pow2(c) * gu * exp(-e0 / (k * T)) / (8 * pi * pow4(f0) * Q);
}

Numeric line::ds_da(Numeric T, Numeric Q) const noexcept {
  return pow2(c) * gu * exp(-e0 / (k * T)) / (8 * pi * pow3(f0) * Q);
}

void band_data::sort(variable v) {
  using enum variable;
  switch (v) {
    case f0:
      std::ranges::sort(lines, {}, &line::f0);
      break;
    case e0:
      std::ranges::sort(lines, {}, &line::e0);
      break;
    case a:
      std::ranges::sort(lines, {}, &line::a);
      break;
    case FINAL:;
  }
}

std::ostream& operator<<(std::ostream& os, const line& x) {
  return os << std::setprecision(std::numeric_limits<Numeric>::digits10 + 1)
            << x.f0 << ' ' << x.a << ' ' << x.e0 << ' ' << x.gu << ' ' << x.gl
            << ' ' << x.z << ' ' << x.ls << ' ' << x.qn.val;
}

std::istream& operator>>(std::istream& is, line& x) {
  is >> double_imanip() >> x.f0 >> x.a >> x.e0 >> x.gu >> x.gl;
  return is >> x.z >> x.ls >> x.qn.val;
}

std::ostream& operator<<(std::ostream& os, const std::vector<line>& x) {
  constexpr std::string_view endl = "\n";
  std::string_view sep = "";
  for (auto& l : x) {
    os << sep << l;
    std::exchange(sep, endl);
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const band_data& x) {
  return os << x.lineshape << ' ' << x.cutoff << ' ' << x.cutoff_value << '\n'
            << x.lines;
}

std::ostream& operator<<(std::ostream& os, const band& x) {
  return os << x.key << '\n' << x.data;
}

std::ostream& operator<<(std::ostream& os, const std::vector<band>& x) {
  constexpr std::string_view endl = "\n";
  std::string_view sep = "";
  for (auto& y : x) {
    os << sep << y;
    std::exchange(sep, endl);
  }
  return os;
}

//! Gets all the lines between (f0-cutoff, f1+cutoff)
std::span<const line> band_data::active_lines(Numeric f0, Numeric f1) const {
  const Numeric c = get_cutoff_frequency();
  auto low = std::ranges::lower_bound(*this, f0 - c, {}, &line::f0);
  auto upp = std::ranges::upper_bound(low, end(), f1 + c, {}, &line::f0);

  return {low, upp};
}

std::ostream& operator<<(std::ostream& os, const line_key& x) {
  return os << "line_key:\n  band: " << x.band << "\n  line: " << x.line
            << "\n  spec: " << x.spec << "\n  var: " << x.var
            << "\n  ls_var: " << x.ls_var << "\n  ls_coeff: " << x.ls_coeff
            << '\n';
}
}  // namespace lbl
