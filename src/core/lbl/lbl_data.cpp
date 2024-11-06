#include "lbl_data.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include "arts_constants.h"
#include "arts_constexpr_math.h"
#include "debug.h"
#include "double_imanip.h"
#include "hitran_species.h"
#include "partfun.h"
#include "quantum_numbers.h"

//! In CPP file
using Constant::c;
using Constant::h;
using Constant::k;
using Constant::pi;
using Math::pow2;
using Math::pow3;
using Math::pow4;
using std::exp;
using std::expm1;

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
    case f0: std::ranges::sort(lines, {}, &line::f0); break;
    case e0: std::ranges::sort(lines, {}, &line::e0); break;
    case a:  std::ranges::sort(lines, {}, &line::a); break;
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
  ARTS_USER_ERROR_IF(not x.qn.val.good(), "Bad quantum numbers in {}", x.qn)

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

std::ostream& operator<<(std::ostream& os, const AbsorptionBands& x) {
  constexpr std::string_view endl = "\n";
  std::string_view sep            = "";
  for (auto& [key, data] : x) {
    os << sep << key << '\n' << data;
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
  auto ptr = absorption_bands.find(type.band);

  ARTS_USER_ERROR_IF(ptr == absorption_bands.end(),
                     "No band with quantum identifier: {}",
                     type.band);

  auto& band = ptr->second;

  ARTS_USER_ERROR_IF(type.line >= band.lines.size(),
                     "Line index out of range: {}"
                     " band has {}"
                     " absorption lines. Band: {}",
                     type.line,
                     band.lines.size(),
                     type.band);
  auto& line = band.lines[type.line];

  if (good_enum(type.ls_var)) {
    auto& line_ls_data = line.ls.single_models;
    ARTS_USER_ERROR_IF(type.spec >= line_ls_data.size(),
                       "Not enough line data for line: {}"
                       " for quantum identifier: {}",
                       line,
                       type.band);

    auto& ls_data = line_ls_data[type.spec].data;
    auto ls_ptr   = std::ranges::find_if(
        ls_data, [var = type.ls_var](auto& x) { return x.first == var; });
    ARTS_USER_ERROR_IF(ls_ptr == ls_data.end(),
                       "No line shape parameter: \"{}"
                       "\" for species: \"{}"
                       "\" in line: {}"
                       " for quantum identifier: {}",
                       type.ls_var,
                       line_ls_data[type.spec].species,
                       line,
                       type.band);

    return ls_ptr->second.X(type.ls_coeff);
  }

  switch (type.var) {
    case LineByLineVariable::f0: return line.f0;
    case LineByLineVariable::e0: return line.e0;
    case LineByLineVariable::a:  return line.a;
  }

  std::unreachable();
}

Numeric& line_key::get_value(AbsorptionBands& b) const {
  return local_get_value(b, *this);
}

const Numeric& line_key::get_value(const AbsorptionBands& b) const {
  return local_get_value(b, *this);
}

Numeric line::hitran_a(const Numeric hitran_s,
                       const SpeciesIsotope& isot,
                       const Numeric T0) const {
  const Numeric Q0 = PartitionFunctions::Q(T0, isot);
  const Numeric Ia = Hitran::isotopologue_ratios()[isot];

  //! Note negative value because expm1 is used as a more accurate form of (1 - exp(x)) for exp(x) close to 1.
  return -8.0 * pi * Q0 * hitran_s /
         (Ia * gu * exp(-e0 / (k * T0)) * expm1(-(h * f0) / (k * T0)) *
          pow2(c / f0));
}

Numeric line::hitran_s(const SpeciesIsotope& isot, const Numeric T0) const {
  return a / hitran_a(1.0, isot, T0);
}

bool band_data::merge(const line& linedata) {
  for (auto& line : lines) {
    if (line.qn == linedata.qn) {
      line = linedata;
      return false;
    }
  }
  lines.push_back(linedata);
  return true;
}

std::unordered_set<SpeciesEnum> species_in_bands(
    const std::unordered_map<QuantumIdentifier, band_data>& bands) {
  std::unordered_set<SpeciesEnum> out;
  for (auto& [key, data] : bands) {
    out.insert(key.Species());

    for (auto& line : data.lines) {
      for (auto spec :
           line.ls.single_models |
               std::views::transform([](auto& x) { return x.species; })) {
        out.insert(spec);
      }
    }
  }
  return out;
}

void keep_hitran_s(std::unordered_map<QuantumIdentifier, band_data>& bands,
                   const std::unordered_map<SpeciesEnum, Numeric>& keep,
                   const Numeric T0) {
  for (auto& [key, data] : bands) {
    const auto ptr = keep.find(key.Species());
    if (ptr != keep.end()) {
      std::erase_if(data.lines, [&key, &T0, &min_s = ptr->second](line& line) {
        return line.hitran_s(key.Isotopologue(), T0) < min_s;
      });
    }
  }
}

std::unordered_map<SpeciesEnum, Numeric> percentile_hitran_s(
    const std::unordered_map<QuantumIdentifier, band_data>& bands,
    const Numeric approx_percentile,
    const Numeric T0) {
  ARTS_USER_ERROR_IF(approx_percentile < 0 or approx_percentile > 100,
                     "Approximate percentile must be between 0 and 100");

  std::unordered_map<SpeciesEnum, std::vector<Numeric>> compute;
  for (auto& [key, data] : bands) {
    for (auto& line : data.lines) {
      compute[key.Species()].push_back(line.hitran_s(key.Isotopologue(), T0));
    }
  }

  std::unordered_map<SpeciesEnum, Numeric> out;
  for (auto& [spec, values] : compute) {
    if (const Size N = values.size(); N != 0) {
      std::ranges::sort(values);
      const Size i =
          static_cast<Size>(static_cast<Numeric>(N) * approx_percentile * 0.01);
      out[spec] = values[std::clamp<Size>(i, 0, N - 1)];
    }
  }

  return out;
}

std::unordered_map<SpeciesEnum, Numeric> percentile_hitran_s(
    const std::unordered_map<QuantumIdentifier, band_data>& bands,
    const std::unordered_map<SpeciesEnum, Numeric>& approx_percentile,
    const Numeric T0) {
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(approx_percentile | std::views::values,
                          [](auto i) { return i < 0 or i > 100; }),
      "Approximate percentile must be between 0 and 100");

  std::unordered_map<SpeciesEnum, std::vector<Numeric>> compute;
  for (auto& [key, data] : bands) {
    if (approx_percentile.contains(key.Species())) {
      for (auto& line : data.lines) {
        compute[key.Species()].push_back(line.hitran_s(key.Isotopologue(), T0));
      }
    }
  }

  std::unordered_map<SpeciesEnum, Numeric> out;
  for (auto& [spec, values] : compute) {
    if (const Size N = values.size(); N != 0) {
      std::ranges::sort(values);
      const Size i = static_cast<Size>(static_cast<Numeric>(N) *
                                       approx_percentile.at(spec) * 0.01);
      out[spec]    = values[std::clamp<Size>(i, 0, N - 1)];
    }
  }

  return out;
}
}  // namespace lbl
