#include "lbl_data.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <type_traits>
#include <utility>

#include "arts_constants.h"
#include "arts_constexpr_math.h"
#include "debug.h"
#include "double_imanip.h"
#include "quantum_numbers.h"

//! In CPP file
using Constant::k;
using Math::pow2;
using Math::pow3;
using Math::pow4;
using std::exp;

namespace lbl {
Numeric line::s(Numeric T, Numeric Q) const {
  return a * gu * exp(-e0 / (k * T)) / (pow3(f0) * Q);
}

Numeric line::ds_dT(Numeric T, Numeric Q, Numeric dQ_dT) const {
  return a * gu * (e0 * Q - k * pow2(T) * dQ_dT) * exp(-e0 / (k * T)) /
         (pow3(f0) * k * pow2(T) * pow2(Q));
}

Numeric line::ds_de0(Numeric T, Numeric Q) const {
  return -a * gu * exp(-e0 / (k * T)) / (pow3(f0) * k * T * Q);
}

Numeric line::ds_df0(Numeric T, Numeric Q) const {
  return -3 * a * gu * exp(-e0 / (k * T)) / (pow4(f0) * Q);
}

Numeric line::ds_da(Numeric T, Numeric Q) const {
  return gu * exp(-e0 / (k * T)) / (pow3(f0) * Q);
}

void band_data::sort(LineByLineVariable v) {
  using enum LineByLineVariable;
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
  }
}

std::ostream& operator<<(std::ostream& os, const line& x) {
  return os << std::setprecision(std::numeric_limits<Numeric>::digits10 + 1)
            << x.f0 << ' ' << x.a << ' ' << x.e0 << ' ' << x.gu << ' ' << x.gl
            << ' ' << x.z << ' ' << x.ls << ' ' << x.qn.val.size() << ' '
            << x.qn.val;
}

std::istream& operator>>(std::istream& is, line& x) {
  Size s{};

  is >> double_imanip() >> x.f0 >> x.a >> x.e0 >> x.gu >> x.gl;
  is >> x.z >> x.ls >> s;
  x.qn.val.reserve(s);
  for (Size i = 0; i < s; i++) is >> x.qn.val.emplace_back();
  ARTS_USER_ERROR_IF(not x.qn.val.good(), "Bad quantum numbers in ", x.qn)

  return is;
}

std::ostream& operator<<(std::ostream& os, const std::vector<line>& x) {
  constexpr std::string_view endl = "\n";
  std::string_view sep            = "";
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
  std::string_view sep            = "";
  for (auto& y : x) {
    os << sep << y;
    std::exchange(sep, endl);
  }
  return os;
}

//! Gets all the lines between (f0-cutoff, f1+cutoff) and the offset from the front
std::pair<Size, std::span<const line>> band_data::active_lines(
    Numeric f0, Numeric f1) const {
  const Numeric c = get_cutoff_frequency();
  auto low        = std::ranges::lower_bound(*this, f0 - c, {}, &line::f0);
  auto upp        = std::ranges::upper_bound(low, end(), f1 + c, {}, &line::f0);

  return {static_cast<Size>(std::distance(begin(), low)), {low, upp}};
}

Rational band_data::max(QuantumNumberType x) const {
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(
          lines, [x](auto& qn) { return not qn.val.has(x); }, &line::qn),
      "Quantum number not found in all lines of the band");

  Rational out{std::numeric_limits<Index>::lowest()};
  for (auto& line : *this) {
    auto& qn = line.qn.val[x];
    out      = std::max(out, std::max(qn.upp(), qn.low()));
  }
  return out;
}

std::ostream& operator<<(std::ostream& os, const line_key& x) {
  return os << "line_key:\n  band: " << x.band << "\n  line: " << x.line
            << "\n  spec: " << x.spec << "\n  var: " << x.var
            << "\n  ls_var: " << x.ls_var << "\n  ls_coeff: " << x.ls_coeff
            << '\n';
}

template <typename T>
auto local_get_value(T& absorption_bands, const line_key& type)
    -> std::conditional_t<std::is_const_v<T>, const Numeric&, Numeric&> {
  auto& band =
      [&type, &absorption_bands]() {
        auto ptr =
            std::ranges::find(absorption_bands, type.band, &lbl::band::key);
        ARTS_USER_ERROR_IF(ptr == absorption_bands.end(),
                           "No band with quantum identifier: ",
                           type.band);
        return ptr;
      }()
          ->data;

  ARTS_USER_ERROR_IF(type.line >= band.lines.size(),
                     "Line index out of range: ",
                     type.line,
                     " band has ",
                     band.lines.size(),
                     " absorption lines. Band: ",
                     type.band);
  auto& line = band.lines[type.line];

  if (good_enum(type.ls_var)) {
    auto& line_ls_data = line.ls.single_models;
    ARTS_USER_ERROR_IF(type.spec >= line_ls_data.size(),
                       "Not enough line data for line: ",
                       line,
                       " for quantum identifier: ",
                       type.band);

    auto& ls_data = line_ls_data[type.spec].data;
    auto ls_ptr   = std::ranges::find_if(
        ls_data, [var = type.ls_var](auto& x) { return x.first == var; });
    ARTS_USER_ERROR_IF(ls_ptr == ls_data.end(),
                       "No line shape parameter: \"",
                       type.ls_var,
                       "\" for species: \"",
                       line_ls_data[type.spec].species,
                       "\" in line: ",
                       line,
                       " for quantum identifier: ",
                       type.band);

    return ls_ptr->second.X(type.ls_coeff);
  }

  switch (type.var) {
    case LineByLineVariable::f0:
      return line.f0;
    case LineByLineVariable::e0:
      return line.e0;
    case LineByLineVariable::a:
      return line.a;
  }

  std::unreachable();
}

Numeric& line_key::get_value(std::vector<lbl::band>& b) const {
  return local_get_value(b, *this);
}

const Numeric& line_key::get_value(const std::vector<lbl::band>& b) const {
  return local_get_value(b, *this);
}
}  // namespace lbl
