#include "lbl_data.h"

#include <arts_constants.h>
#include <arts_constexpr_math.h>
#include <debug.h>
#include <double_imanip.h>
#include <hitran_species.h>
#include <partfun.h>
#include <quantum.h>
#include <rational.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <format>
#include <limits>
#include <print>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>

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

std::istream& operator>>(std::istream& is, line& x) try {
  is >> double_imanip() >> x.f0 >> x.a >> x.e0 >> x.gu >> x.gl;
  return is >> x.z >> x.ls >> x.qn;
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error reading line data:\n{}\n{:IO}", e.what(), x));
}

//! Gets all the lines between (f0-cutoff, f1+cutoff) and the offset from the front
std::pair<Size, std::span<const line>> band_data::active_lines(
    Numeric f0, Numeric f1) const {
  const Numeric c = get_cutoff_frequency();
  auto low        = std::ranges::lower_bound(*this, f0 - c, {}, &line::f0);
  auto upp        = std::ranges::upper_bound(low, end(), f1 + c, {}, &line::f0);

  return {static_cast<Size>(std::distance(begin(), low)), {low, upp}};
}

Rational band_data::max(QuantumNumberType x) const try {
  ARTS_USER_ERROR_IF(lines.empty(), "Must have lines")

  Rational out{maxr(lines.front().qn.at(x).upper.get<Rational>(),
                    lines.front().qn.at(x).lower.get<Rational>())};
  for (auto& line : lines | stdv::drop(1)) {
    auto& qn = line.qn.at(x);
    out = maxr(out, maxr(qn.upper.get<Rational>(), qn.lower.get<Rational>()));
  }
  return out;
} catch (std::out_of_range&) {
  throw std::runtime_error(std::format(
      "Cannot find QuantumNumberType \"{}\" for at least on line", x));
} catch (std::exception&) {
  throw;
}

namespace {
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
}  // namespace

Numeric& line_key::get_value(AbsorptionBands& b) const {
  return local_get_value(b, *this);
}

const Numeric& line_key::get_value(const AbsorptionBands& b) const {
  return local_get_value(b, *this);
}

[[nodiscard]] Numeric line::compute_a(const Numeric s,
                                      const SpeciesIsotope& isot,
                                      const Numeric T0) const {
  const Numeric Q0 = PartitionFunctions::Q(T0, isot);

  //! Note negative value because expm1 is used as a more accurate form of (1 - exp(x)) for exp(x) close to 1.
  return -8.0 * pi * Q0 * s /
         (gu * exp(-e0 / (k * T0)) * expm1(-(h * f0) / (k * T0)) *
          pow2(c / f0));
}

Numeric line::hitran_a(const Numeric hitran_s,
                       const SpeciesIsotope& isot,
                       const Numeric T0) const {
  const Numeric Ia = Hitran::isotopologue_ratios()[isot];
  return compute_a(hitran_s / Ia, isot, T0);
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
    out.insert(key.isot.spec);

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
    const auto ptr = keep.find(key.isot.spec);
    if (ptr != keep.end()) {
      std::erase_if(data.lines, [&key, &T0, &min_s = ptr->second](line& line) {
        return line.hitran_s(key.isot, T0) < min_s;
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
      compute[key.isot.spec].push_back(line.hitran_s(key.isot, T0));
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
    if (approx_percentile.contains(key.isot.spec)) {
      for (auto& line : data.lines) {
        compute[key.isot.spec].push_back(line.hitran_s(key.isot, T0));
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

Size count_lines(
    const std::unordered_map<QuantumIdentifier, lbl::band_data>& bands) {
  return std::transform_reduce(
      bands.begin(), bands.end(), Size{0}, std::plus{}, [](const auto& x) {
        return x.second.size();
      });
}
}  // namespace lbl

std::string std::formatter<lbl::line>::to_string(const lbl::line& v) const {
  if (tags.io) {
    return tags.vformat(v.f0,
                        ' ',
                        v.a,
                        ' ',
                        v.e0,
                        ' ',
                        v.gu,
                        ' ',
                        v.gl,
                        ' ',
                        v.z,
                        ' ',
                        v.ls,
                        ' ',
                        v.qn);
  }

  const std::string_view sep = tags.sep();
  std::string out =
      tags.vformat(v.f0, sep, v.a, sep, v.e0, sep, v.gu, sep, v.gl);
  if (not tags.short_str) out += tags.vformat(sep, v.z, sep, v.ls, sep, v.qn);

  return tags.bracket ? ("[" + out + "]") : out;
}

std::string std::formatter<lbl::band_data>::to_string(
    const lbl::band_data& v) const {
  if (tags.io) {
    return tags.vformat(v.lineshape, ' ', v.cutoff, ' ', v.cutoff_value);
  }

  const auto sep = tags.sep();

  std::string out =
      tags.vformat(v.lineshape, sep, v.cutoff, sep, v.cutoff_value);
  if (not tags.short_str) out += tags.vformat(sep, v.lines);

  return out;
}

void xml_io_stream<AbsorptionBand>::write(std::ostream& os,
                                          const AbsorptionBand& data,
                                          bofstream* pbofs,
                                          std::string_view name) {
  ARTS_USER_ERROR_IF(pbofs not_eq nullptr, "No binary data")

  XMLTag open_tag;
  open_tag.name=("AbsorptionBandData");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("lineshape", String{toString(data.lineshape)});
  open_tag.add_attribute("cutoff_type", String{toString(data.cutoff)});
  open_tag.add_attribute("cutoff_value", data.cutoff_value);
  open_tag.add_attribute("nelem", static_cast<Index>(data.lines.size()));
  open_tag.write_to_stream(os);
  std::println(os);

  for (auto& line : data) {
    std::println(os, "{:IO}", line);
  }

  XMLTag close_tag;
  close_tag.name=("/AbsorptionBandData");
  close_tag.write_to_stream(os);
  std::println(os);
}

void xml_io_stream<AbsorptionBand>::read(std::istream& is,
                                         AbsorptionBand& data,
                                         bifstream* pbifs) try {
  ARTS_USER_ERROR_IF(pbifs not_eq nullptr, "No binary data")

  String tag;
  Index nelem;

  XMLTag open_tag;
  open_tag.read_from_stream(is);
  open_tag.check_name("AbsorptionBandData");

  open_tag.get_attribute_value("lineshape", tag);
  data.lineshape = to<LineByLineLineshape>(tag);

  open_tag.get_attribute_value("cutoff_type", tag);
  data.cutoff = to<LineByLineCutoffType>(tag);

  open_tag.get_attribute_value("cutoff_value", data.cutoff_value);

  open_tag.get_attribute_value("nelem", nelem);
  data.lines.resize(0);
  data.lines.reserve(nelem);

  for (Index j = 0; j < nelem; j++) {
    is >> data.lines.emplace_back();
  }

  XMLTag close_tag;
  close_tag.read_from_stream(is);
  close_tag.check_name("/AbsorptionBandData");
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error reading AbsorptionBandData:\n{}", e.what()));
}

void xml_io_stream<LblLineKey>::write(std::ostream& os,
                                      const LblLineKey& x,
                                      bofstream* pbofs,
                                      std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.band, pbofs);
  xml_write_to_stream(os, x.line, pbofs);
  xml_write_to_stream(os, x.spec, pbofs);
  xml_write_to_stream(os, x.ls_var, pbofs);
  xml_write_to_stream(os, x.ls_coeff, pbofs);
  xml_write_to_stream(os, x.var, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<LblLineKey>::read(std::istream& is,
                                     LblLineKey& x,
                                     bifstream* pbifs) try {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.band, pbifs);
  xml_read_from_stream(is, x.line, pbifs);
  xml_read_from_stream(is, x.spec, pbifs);
  xml_read_from_stream(is, x.ls_var, pbifs);
  xml_read_from_stream(is, x.ls_coeff, pbifs);
  xml_read_from_stream(is, x.var, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error reading LblLineKey:\n{}", e.what()));
}
