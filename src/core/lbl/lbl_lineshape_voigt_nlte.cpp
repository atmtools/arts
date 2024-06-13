#include "lbl_lineshape_voigt_nlte.h"

#include <jacobian.h>
#include <partfun.h>
#include <physics_funcs.h>
#include <sorting.h>

#include <Faddeeva/Faddeeva.hh>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>

#include "arts_constants.h"
#include "arts_constexpr_math.h"
#include "atm.h"
#include "debug.h"
#include "lbl_data.h"
#include "lbl_zeeman.h"
#include "quantum_numbers.h"
#include "rtepack.h"
#include "species.h"

namespace lbl::voigt::nlte {
std::pair<Numeric, Numeric> line_strength_calc(const Numeric inv_gd,
                                               const QuantumIdentifier& qid,
                                               const line& line,
                                               const AtmPoint& atm) {
  using Constant::h, Constant::inv_sqrt_pi;
  using Math::pow2, Math::pow3;

  const Numeric ru = atm[qid.UpperLevel()];
  const Numeric rl = atm[qid.LowerLevel()];

  const Numeric k             = line.nlte_k(ru, rl);
  const Numeric e             = line.nlte_e(ru);
  const Numeric f             = line.f0;
  static constexpr Numeric kB = Constant::k;
  const Numeric T             = atm.temperature;
  const Numeric inv_b         = -std::expm1(h * f / (kB * T)) / (f * f * f);

  const Numeric r = atm[qid.Isotopologue()];
  const Numeric x = atm[qid.Isotopologue().spec];

  //! Missing factor is c^2 f / 8pi
  return {inv_sqrt_pi * inv_gd * r * x * k,
          inv_sqrt_pi * inv_gd * r * x * (e * inv_b - k)};
}

std::pair<Numeric, Numeric> dline_strength_calc_dT(const Numeric inv_gd,
                                                   const Numeric f0,
                                                   const QuantumIdentifier& qid,
                                                   const line& line,
                                                   const AtmPoint& atm) {
  using Constant::h, Constant::inv_sqrt_pi;
  using Math::pow2, Math::pow3;

  const Numeric ru = atm[qid.UpperLevel()];
  const Numeric rl = atm[qid.LowerLevel()];

  const Numeric k             = line.nlte_k(ru, rl);
  const Numeric e             = line.nlte_e(ru);
  const Numeric f             = line.f0;
  static constexpr Numeric kB = Constant::k;
  const Numeric T             = atm.temperature;
  const Numeric inv_b         = -std::expm1(h * f / (kB * T)) / (f * f * f);
  const Numeric dinv_b        = -h * exp(f * h / (T * k)) / (T * T * f * f * k);

  const Numeric dD0 = line.ls.dD0_dT(atm);

  const Numeric dinv_gd = inv_gd * (2 * T * dD0 + f0) / (2 * T * f0);

  const Numeric r = atm[qid.Isotopologue()];
  const Numeric x = atm[qid.Isotopologue().spec];

  //! Missing factor is c^2 f / 8pi
  return {inv_sqrt_pi * dinv_gd * r * x * k,
          inv_sqrt_pi * inv_gd * r * x * e * dinv_b +
              inv_sqrt_pi * dinv_gd * r * x * (e * inv_b - k)};
}

std::pair<Numeric, Numeric> dline_strength_calc_dT(const Numeric inv_gd,
                                                   const Numeric f0,
                                                   const QuantumIdentifier& qid,
                                                   const line& line,
                                                   const AtmPoint& atm,
                                                   const Size ispec) {
  using Constant::h, Constant::inv_sqrt_pi;
  using Math::pow2, Math::pow3;

  const Numeric ru = atm[qid.UpperLevel()];
  const Numeric rl = atm[qid.LowerLevel()];

  const Numeric k             = line.nlte_k(ru, rl);
  const Numeric e             = line.nlte_e(ru);
  const Numeric f             = line.f0;
  static constexpr Numeric kB = Constant::k;
  const Numeric T             = atm.temperature;
  const Numeric inv_b         = -std::expm1(h * f / (kB * T)) / (f * f * f);
  const Numeric dinv_b        = -h * exp(f * h / (T * k)) / (T * T * f * f * k);

  const auto& ls    = line.ls.single_models[ispec];
  const Numeric T0  = line.ls.T0;
  const Numeric P   = atm.pressure;
  const Numeric dD0 = ls.dD0_dT(T0, T, P);

  const Numeric dinv_gd = inv_gd * (2 * T * dD0 + f0) / (2 * T * f0);

  const Numeric r = atm[qid.Isotopologue()];
  const Numeric x = atm[qid.Isotopologue().spec];

  //! Missing factor is c^2 f / 8pi
  return {inv_sqrt_pi * dinv_gd * r * x * k,
          inv_sqrt_pi * inv_gd * r * x * e * dinv_b +
              inv_sqrt_pi * dinv_gd * r * x * (e * inv_b - k)};
}

Numeric line_center_calc(const line& line, const AtmPoint& atm) {
  return line.f0 + line.ls.D0(atm);
}

Numeric dline_center_calc_dT(const line& line, const AtmPoint& atm) {
  return line.ls.dD0_dT(atm);
}

Numeric dline_center_calc_dVMR(const line& line,
                               const SpeciesEnum spec,
                               const AtmPoint& atm) {
  return line.ls.dD0_dVMR(atm, spec);
}

Numeric line_center_calc(const line& line, const AtmPoint& atm, Size ispec) {
  const auto& ls = line.ls.single_models[ispec];
  return line.f0 + ls.D0(line.ls.T0, atm.temperature, atm.pressure);
}

Numeric dline_center_calc_dT(const line& line,
                             const AtmPoint& atm,
                             Size ispec) {
  const auto& ls = line.ls.single_models[ispec];
  return ls.dD0_dT(line.ls.T0, atm.temperature, atm.pressure);
}

Numeric dline_center_calc_dVMR(const line& line,
                               const SpeciesEnum spec,
                               const AtmPoint& atm,
                               Size ispec) {
  const auto& ls = line.ls.single_models[ispec];
  return ls.species == spec ? ls.D0(line.ls.T0, atm.temperature, atm.pressure)
                            : 0;
}

Numeric scaled_gd(const Numeric T, const Numeric mass, const Numeric f0) {
  constexpr auto c = Constant::doppler_broadening_const_squared;
  return std::sqrt(c * T / mass) * f0;
}

//! Should only live in CC-file since it holds references
struct single_shape_builder {
  const QuantumIdentifier& qid;
  const line& ln;
  const AtmPoint& atm;
  Numeric f0;
  Numeric scaled_gd_part;
  Numeric G0;
  Size ispec{std::numeric_limits<Size>::max()};

  single_shape_builder(const QuantumIdentifier& id,
                       const line& l,
                       const AtmPoint& a)
      : qid(id),
        ln(l),
        atm(a),
        f0(line_center_calc(ln, atm)),
        scaled_gd_part(std::sqrt(Constant::doppler_broadening_const_squared *
                                 atm.temperature / id.Isotopologue().mass)),
        G0(ln.ls.G0(atm)) {}

  single_shape_builder(const QuantumIdentifier& id,
                       const line& l,
                       const AtmPoint& a,
                       const Size is)
      : qid(id),
        ln(l),
        atm(a),
        f0(line_center_calc(ln, atm, is)),
        scaled_gd_part(std::sqrt(Constant::doppler_broadening_const_squared *
                                 atm.temperature / id.Isotopologue().mass)),
        G0(ln.ls.single_models[is].G0(ln.ls.T0, atm.temperature, atm.pressure)),
        ispec(is) {}

  [[nodiscard]] single_shape as_zeeman(const Numeric H,
                                       const zeeman::pol pol,
                                       const Size iz) const {
    single_shape s;
    s.f0                    = f0 + H * ln.z.Splitting(ln.qn.val, pol, iz);
    s.inv_gd                = 1.0 / (scaled_gd_part * f0);
    s.z_imag                = G0 * s.inv_gd;
    const auto [k, e_ratio] = line_strength_calc(s.inv_gd, qid, ln, atm);
    s.k                     = ln.z.Strength(ln.qn.val, pol, iz) * k;
    s.e_ratio               = ln.z.Strength(ln.qn.val, pol, iz) * e_ratio;
    return s;
  }

  operator single_shape() const {
    single_shape s;
    s.f0                    = f0;
    s.inv_gd                = 1.0 / (scaled_gd_part * f0);
    s.z_imag                = G0 * s.inv_gd;
    const auto [k, e_ratio] = line_strength_calc(s.inv_gd, qid, ln, atm);
    s.k                     = k;
    s.e_ratio               = e_ratio;
    return s;
  }
};

single_shape::single_shape(const QuantumIdentifier& qid,
                           const line& line,
                           const AtmPoint& atm,
                           const zeeman::pol pol,
                           const Index iz)
    : f0(line_center_calc(line, atm) +
         std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]) *
             line.z.Splitting(line.qn.val, pol, iz)),
      inv_gd(1.0 / scaled_gd(atm.temperature, qid.Isotopologue().mass, f0)),
      z_imag(line.ls.G0(atm) * inv_gd) {
  const auto [kp, e_ratiop] = line_strength_calc(inv_gd, qid, line, atm);
  k                         = line.z.Strength(line.qn.val, pol, iz) * kp;
  e_ratio                   = line.z.Strength(line.qn.val, pol, iz) * e_ratiop;
}

single_shape::single_shape(const QuantumIdentifier& qid,
                           const line& line,
                           const AtmPoint& atm,
                           const zeeman::pol pol,
                           const Index iz,
                           const Size ispec)
    : f0(line_center_calc(line, atm, ispec) +
         std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]) *
             line.z.Splitting(line.qn.val, pol, iz)),
      inv_gd(1.0 / scaled_gd(atm.temperature, qid.Isotopologue().mass, f0)),
      z_imag(line.ls.single_models[ispec].G0(
                 line.ls.T0, atm.temperature, atm.pressure) *
             inv_gd) {
  const auto [kp, e_ratiop] = line_strength_calc(inv_gd, qid, line, atm);
  k                         = line.z.Strength(line.qn.val, pol, iz) * kp;
  e_ratio                   = line.z.Strength(line.qn.val, pol, iz) * e_ratiop;
}

Complex single_shape::F(const Complex z_) { return Faddeeva::w(z_); }

Complex single_shape::F(const Numeric f) const { return F(z(f)); }

std::pair<Complex, Complex> single_shape::operator()(const Numeric f) const {
  const auto F_ = F(f);
  return {k * F_, e_ratio * F_};
}

Complex single_shape::dF(const Numeric f) const {
  const Complex z_ = z(f);
  return dF(z_, F(z_));
}

Complex single_shape::dF(const Complex z_, const Complex F_) {
  /*! FIXME: We should use a proper algorithm here.  This produces
   *         no errors in tests, but the actual derivative form is
   *         analytically known.  Its numerical instability, however,
   *         makes it completely useless.
   *
   *         The analytical form is:
   *
   *           dF = -2 * z * F(z) + 2 * i / sqrt(pi)
   *
   * Tests show that for y < 1e7 it works until x > 1e7, but for
   * y > 1e7, it always fails.  This is about the analytical form
   * above using the latest version of the MIT Faddeeva package.
  */
  const Complex dz{std::max(1e-4 * z_.real(), 1e-4),
                   std::max(1e-4 * z_.imag(), 1e-4)};
  const Complex F_2 = Faddeeva::w(z_ + dz);
  return (F_2 - F_) / dz;
}

single_shape::zFdF::zFdF(const Complex z_)
    : z{z_}, F{single_shape::F(z_)}, dF{single_shape::dF(z_, F)} {}

single_shape::zFdF single_shape::all(const Numeric f) const { return z(f); }
std::pair<Complex, Complex> single_shape::dru(const Numeric dk_dru,
                                              const Numeric de_ratio_dru,
                                              const Numeric f) const {
  const auto F_ = F(f);
  return {dk_dru * F_, de_ratio_dru * F_};
}

std::pair<Complex, Complex> single_shape::drl(const Numeric dk_drl,
                                              const Numeric de_ratio_drl,
                                              const Numeric f) const {
  const auto F_ = F(f);
  return {dk_drl * F_, de_ratio_drl * F_};
}

std::pair<Complex, Complex> single_shape::df(const Numeric f) const {
  const auto dF_ = inv_gd * dF(f);
  return {k * dF_, e_ratio * dF_};
}

std::pair<Complex, Complex> single_shape::dH(const Complex dz_dH,
                                             const Numeric f) const {
  const auto dF_ = dz_dH * dF(f);
  return {k * dF_, e_ratio * dF_};
}

std::pair<Complex, Complex> single_shape::dT(const Numeric dk_dT,
                                             const Numeric de_ratio_dT,
                                             const Complex dz_dT,
                                             const Numeric dz_dT_fac,
                                             const Numeric f) const {
  const auto [z_, F_, dF_] = all(f);
  return {dk_dT * F_ + k * (dz_dT + dz_dT_fac * z_) * dF_,
          de_ratio_dT * F_ + e_ratio * (dz_dT + dz_dT_fac * z_) * dF_};
}

constexpr std::pair<Index, Index> find_offset_and_count_of_frequency_range(
    const std::span<const single_shape> lines, Numeric f, Numeric cutoff) {
  if (cutoff < std::numeric_limits<Numeric>::infinity()) {
    auto low =
        std::ranges::lower_bound(lines, f - cutoff, {}, &single_shape::f0);
    auto upp =
        std::ranges::upper_bound(lines, f + cutoff, {}, &single_shape::f0);

    return {std::distance(lines.begin(), low), std::distance(low, upp)};
  }

  return {0, lines.size()};
}

namespace detail {
constexpr auto frequency_span(const auto& list,
                              const Size start,
                              const Size count) {
  std::span out{list};
  return out.subspan(start, count);
}
}  // namespace detail

template <typename... Ts>
constexpr auto frequency_spans(const Numeric cutoff,
                               const Numeric f,
                               const std::span<const single_shape>& lines,
                               const Ts&... lists) {
  ARTS_ASSERT(lines.size() == (static_cast<Size>(lists.size()) and ...))

  const auto [start, count] =
      find_offset_and_count_of_frequency_range(lines, f, cutoff);
  return std::tuple{detail::frequency_span(lines, start, count),
                    detail::frequency_span(lists, start, count)...};
}

Size count_lines(const band_data& bnd, const zeeman::pol type) {
  return std::transform_reduce(
      bnd.begin(), bnd.end(), Index{}, std::plus<>{}, [type](auto& line) {
        const Index factor =
            line.ls.one_by_one ? line.ls.single_models.size() : 1;
        return factor * line.z.size(line.qn.val, type);
      });
}

void zeeman_push_back(std::vector<single_shape>& lines,
                      std::vector<line_pos>& pos,
                      const single_shape_builder& s,
                      const line& line,
                      const AtmPoint& atm,
                      const zeeman::pol pol,
                      const Size ispec,
                      const Size iline) {
  if (pol == zeeman::pol::no) {
    lines.emplace_back(s);
    pos.emplace_back(line_pos{.line = iline, .spec = ispec});
  } else {
    const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
    const auto nz   = static_cast<Size>(line.z.size(line.qn.val, pol));
    for (Size iz = 0; iz < nz; iz++) {
      lines.emplace_back(s.as_zeeman(H, pol, iz));
      pos.emplace_back(line_pos{.line = iline, .spec = ispec, .iz = iz});
    }
  }
}

void lines_push_back(std::vector<single_shape>& lines,
                     std::vector<line_pos>& pos,
                     const QuantumIdentifier& qid,
                     const line& line,
                     const AtmPoint& atm,
                     const zeeman::pol pol,
                     const Size iline) {
  if (line.ls.one_by_one) {
    for (Size i = 0; i < line.ls.single_models.size(); ++i) {
      if ((line.z.on and pol != zeeman::pol::no) or
          (not line.z.on and pol == zeeman::pol::no)) {
        zeeman_push_back(lines,
                         pos,
                         single_shape_builder{qid, line, atm, i},
                         line,
                         atm,
                         pol,
                         i,
                         iline);
      }
    }
  } else {
    if ((line.z.on and pol != zeeman::pol::no) or
        (not line.z.on and pol == zeeman::pol::no)) {
      zeeman_push_back(lines,
                       pos,
                       single_shape_builder{qid, line, atm},
                       line,
                       atm,
                       pol,
                       std::numeric_limits<Size>::max(),
                       iline);
    }
  }
}

void band_shape_helper(std::vector<single_shape>& lines,
                       std::vector<line_pos>& pos,
                       const QuantumIdentifier& qid,
                       const band_data& bnd,
                       const AtmPoint& atm,
                       const Numeric fmin,
                       const Numeric fmax,
                       const zeeman::pol pol) {
  lines.resize(0);
  pos.resize(0);

  lines.reserve(count_lines(bnd, pol));
  pos.reserve(lines.capacity());

  using enum LineByLineCutoffType;
  switch (bnd.cutoff) {
    case None:
      for (Size iline = 0; iline < bnd.size(); iline++) {
        lines_push_back(lines, pos, qid, bnd.lines[iline], atm, pol, iline);
      }
      break;
    case ByLine: {
      auto [iline, active_lines] = bnd.active_lines(fmin, fmax);
      for (auto& line : active_lines) {
        lines_push_back(lines, pos, qid, line, atm, pol, iline++);
      }
    } break;
  }

  bubble_sort_by(
      [&](const Size l1, const Size l2) { return lines[l1].f0 > lines[l2].f0; },
      lines,
      pos);
}

band_shape::band_shape(std::vector<single_shape>&& ls, const Numeric cut)
    : lines(std::move(ls)), cutoff(cut) {}

constexpr static auto add_pair = [](auto&& lhs,
                                    auto&& rhs) -> std::pair<Complex, Complex> {
  return {lhs.first + rhs.first, lhs.second + rhs.second};
};

constexpr std::pair<Complex, Complex> rem_pair(auto&& lhs, auto&& rhs) {
  return {lhs.first - rhs.first, lhs.second - rhs.second};
};

std::pair<Complex, Complex> band_shape::operator()(const Numeric f) const {
  return std::transform_reduce(lines.begin(),
                               lines.end(),
                               std::pair<Complex, Complex>{},
                               add_pair,
                               [f](auto& ls) { return ls(f); });
}

std::pair<Complex, Complex> band_shape::df(const Numeric f) const {
  return std::transform_reduce(lines.begin(),
                               lines.end(),
                               std::pair<Complex, Complex>{},
                               add_pair,
                               [f](auto& ls) { return ls.df(f); });
}

std::pair<Complex, Complex> band_shape::dH(
    const ExhaustiveConstComplexVectorView& dz_dH, const Numeric f) const {
  ARTS_ASSERT(static_cast<Size>(dz_dH.size()) == lines.size())

  return std::transform_reduce(lines.begin(),
                               lines.end(),
                               dz_dH.begin(),
                               std::pair<Complex, Complex>{},
                               add_pair,
                               [f](auto& ls, auto& d) { return ls.dH(d, f); });
}

std::pair<Complex, Complex> band_shape::dT(
    const ExhaustiveConstVectorView& dk_dT,
    const ExhaustiveConstVectorView& de_ratio_dT,
    const ExhaustiveConstComplexVectorView& dz_dT,
    const ExhaustiveConstVectorView& dz_dT_fac,
    const Numeric f) const {
  ARTS_ASSERT(dk_dT.size() == dz_dT.size())
  ARTS_ASSERT(static_cast<Size>(dk_dT.size()) == lines.size())

  std::pair<Complex, Complex> out{};  //! Fixme, use zip in C++ 23...

  for (Size i = 0; i < lines.size(); ++i) {
    out = add_pair(
        out, lines[i].dT(dk_dT[i], de_ratio_dT[i], dz_dT[i], dz_dT_fac[i], f));
  }

  return out;
}

std::pair<Complex, Complex> band_shape::operator()(const CutViewConst& cut,
                                                   const Numeric f) const {
  const auto [s, cs] = frequency_spans(cutoff, f, lines, cut);
  return std::transform_reduce(
      s.begin(),
      s.end(),
      cs.begin(),
      std::pair<Complex, Complex>{},
      add_pair,
      [f](auto& ls, auto& c) { return rem_pair(ls(f), c); });
}

void band_shape::operator()(CutView cut) const {
  std::transform(
      lines.begin(),
      lines.end(),
      cut.begin(),
      [cutoff_freq = cutoff](auto& ls) { return ls(ls.f0 + cutoff_freq); });
}

std::pair<Complex, Complex> band_shape::df(const CutViewConst& cut,
                                           const Numeric f) const {
  const auto [s, cs] = frequency_spans(cutoff, f, lines, cut);
  return std::transform_reduce(
      s.begin(),
      s.end(),
      cs.begin(),
      std::pair<Complex, Complex>{},
      add_pair,
      [f](auto& ls, auto& c) { return rem_pair(ls.df(f), c); });
}

void band_shape::df(CutView cut) const {
  std::transform(
      lines.begin(),
      lines.end(),
      cut.begin(),
      [cutoff_freq = cutoff](auto& ls) { return ls.df(ls.f0 + cutoff_freq); });
}

std::pair<Complex, Complex> band_shape::dH(
    const CutViewConst& cut,
    const ExhaustiveConstComplexVectorView& dz_dH,
    const Numeric f) const {
  ARTS_ASSERT(static_cast<Size>(dz_dH.size()) == lines.size())

  const auto [s, cs, dH] = frequency_spans(cutoff, f, lines, cut, dz_dH);

  std::pair<Complex, Complex> out{};  //! Fixme, use zip in C++ 23...
  for (Size i = 0; i < s.size(); ++i) {
    out = add_pair(out, rem_pair(s[i].dH(dH[i], f), cs[i]));
  }

  return out;
}

void band_shape::dH(CutView cut,
                    const ExhaustiveConstComplexVectorView& df0_dH) const {
  ARTS_ASSERT(static_cast<Size>(df0_dH.size()) == lines.size())

  std::transform(lines.begin(),
                 lines.end(),
                 df0_dH.begin(),
                 cut.begin(),
                 [cutoff_freq = cutoff](auto& ls, auto& d) {
                   return ls.dH(d, ls.f0 + cutoff_freq);
                 });
}

std::pair<Complex, Complex> band_shape::dT(
    const CutViewConst& cut,
    const ExhaustiveConstVectorView& dk_dT,
    const ExhaustiveConstVectorView& de_ratio_dT,
    const ExhaustiveConstComplexVectorView& dz_dT,
    const ExhaustiveConstVectorView& dz_dT_fac,
    const Numeric f) const {
  ARTS_ASSERT(dk_dT.size() == dz_dT.size())
  ARTS_ASSERT(static_cast<Size>(dk_dT.size()) == lines.size())

  std::pair<Complex, Complex> out{};  //! Fixme, use zip in C++ 23...

  const auto [s, cs, dk, de, dz, dzf] = frequency_spans(
      cutoff, f, lines, cut, dk_dT, de_ratio_dT, dz_dT, dz_dT_fac);

  for (Size i = 0; i < s.size(); ++i) {
    out =
        add_pair(out, rem_pair(s[i].dT(dk[i], de[i], dz[i], dzf[i], f), cs[i]));
  }

  return out;
}

void band_shape::dT(CutView cut,
                    const ExhaustiveConstVectorView& dk_dT,
                    const ExhaustiveConstVectorView& de_ratio_dT,
                    const ExhaustiveConstComplexVectorView& dz_dT,
                    const ExhaustiveConstVectorView& dz_dT_fac) const {
  ARTS_ASSERT(dk_dT.size() == dz_dT.size())
  ARTS_ASSERT(static_cast<Size>(dk_dT.size()) == lines.size())

  for (Size i = 0; i < lines.size(); ++i) {
    cut[i] = lines[i].dT(
        dk_dT[i], de_ratio_dT[i], dz_dT[i], dz_dT_fac[i], lines[i].f0 + cutoff);
  }
}

void ComputeData::update_zeeman(const Vector2& los,
                                const Vector3& mag,
                                const zeeman::pol pol) {
  npm = zeeman::norm_view(pol, mag, los);
  if (pol != zeeman::pol::no) {
    dnpm_du = zeeman::dnorm_view_du(pol, mag, los);
    dnpm_dv = zeeman::dnorm_view_dv(pol, mag, los);
    dnpm_dw = zeeman::dnorm_view_dw(pol, mag, los);
  }
}

ComputeData::ComputeData(const ExhaustiveConstVectorView& f_grid,
                         const AtmPoint& atm,
                         const Vector2& los,
                         const zeeman::pol pol)
    : scl(f_grid.size()),
      dscl(f_grid.size()),
      shape(f_grid.size()),
      dshape(f_grid.size()) {
  std::transform(f_grid.begin(),
                 f_grid.end(),
                 scl.begin(),
                 [N = number_density(atm.pressure, atm.temperature)](auto f) {
                   constexpr Numeric c =
                       Constant::c * Constant::c / (8 * Constant::pi);
                   return N * f * c;  // Lacking term???
                 });

  update_zeeman(los, atm.mag, pol);
}

//! Sizes cut, dcut, dz, ds; sets shape
void ComputeData::core_calc(const band_shape& shp,
                            const band_data& bnd,
                            const ExhaustiveConstVectorView& f_grid) {
  cut.resize(shp.size());
  dz.resize(shp.size());
  dz_fac.resize(shp.size());
  dk.resize(shp.size());
  de_ratio.resize(shp.size());
  dcut.resize(shp.size());

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp(cut);
    std::transform(
        f_grid.begin(), f_grid.end(), shape.begin(), [this, &shp](Numeric f) {
          return shp(cut, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), shape.begin(), [&shp](Numeric f) {
          return shp(f);
        });
  }
}

//! Sets dshape and dscl and ds and dz
void ComputeData::dt_core_calc(const QuantumIdentifier& qid,
                               const band_shape& shp,
                               const band_data& bnd,
                               const ExhaustiveConstVectorView& f_grid,
                               const AtmPoint& atm,
                               const zeeman::pol pol) {
  std::transform(
      f_grid.begin(),
      f_grid.end(),
      dscl.begin(),
      [dN = dnumber_density_dt(atm.pressure, atm.temperature)](auto f) {
        constexpr Numeric c = Constant::c * Constant::c / (8 * Constant::pi);
        return f * dN * c;
      });

  const Numeric T = atm.temperature;
  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    const auto& lshp = shp.lines[i];

    const Numeric& inv_gd = lshp.inv_gd;
    const Numeric& f0     = lshp.f0;

    if (pos[i].spec == std::numeric_limits<Size>::max()) {
      dz_fac[i] = (-2 * T * line.ls.dD0_dT(atm) - f0) / (2 * T * f0);

      const auto [dkp, dep] =
          dline_strength_calc_dT(inv_gd, f0, qid, line, atm);
      dk[i]       = line.z.Strength(line.qn.val, pol, pos[i].iz) * dkp;
      de_ratio[i] = line.z.Strength(line.qn.val, pol, pos[i].iz) * dep;

      dz[i] = inv_gd *
              Complex{-dline_center_calc_dT(line, atm), line.ls.dG0_dT(atm)};
    } else {
      const auto& ls = line.ls.single_models[pos[i].spec];

      dz_fac[i] =
          (-2 * T * ls.dD0_dT(line.ls.T0, T, atm.pressure) - f0) / (2 * T * f0);

      const auto [dkp, dep] =
          dline_strength_calc_dT(inv_gd, f0, qid, line, atm, pos[i].spec);
      dk[i]       = line.z.Strength(line.qn.val, pol, pos[i].iz) * dkp;
      de_ratio[i] = line.z.Strength(line.qn.val, pol, pos[i].iz) * dep;

      dz[i] = inv_gd * Complex{-ls.dD0_dT(line.ls.T0, T, atm.pressure),
                               ls.dG0_dT(line.ls.T0, T, atm.pressure)};
    }
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dT(dcut, dk, de_ratio, dz, dz_fac);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dT(dcut, dk, de_ratio, dz, dz_fac, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dT(dk, de_ratio, dz, dz_fac, f);
        });
  }
}

//! Sets dshape and dscl
void ComputeData::df_core_calc(const band_shape& shp,
                               const band_data& bnd,
                               const ExhaustiveConstVectorView& f_grid,
                               const AtmPoint& atm) {
  std::transform(f_grid.begin(),
                 f_grid.end(),
                 dscl.begin(),
                 [N = number_density(atm.pressure, atm.temperature),
                  T = atm.temperature](auto f) {
                   constexpr Numeric c =
                       Constant::c * Constant::c / (8 * Constant::pi);
                   const Numeric r = (Constant::h * f) / (Constant::k * T);
                   return N * (r * std::exp(-r) - std::expm1(-r)) * c;
                 });

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.df(dcut);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.df(dcut, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [&shp](Numeric f) {
          return shp.df(f);
        });
  }
}

//! Sets dshape and dz
void ComputeData::dmag_u_core_calc(const band_shape& shp,
                                   const band_data& bnd,
                                   const ExhaustiveConstVectorView& f_grid,
                                   const AtmPoint& atm,
                                   const zeeman::pol pol) {
  const Numeric H         = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  const Numeric dH_dmag_u = atm.mag[0] / H;

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    dz[i]            = -shp.lines[pos[i].line].inv_gd * dH_dmag_u *
            line.z.Splitting(line.qn.val, pol, pos[i].iz);
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dH(dcut, dz);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dH(dcut, dz, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dH(dz, f);
        });
  }
}

//! Sets dshape and dz
void ComputeData::dmag_v_core_calc(const band_shape& shp,
                                   const band_data& bnd,
                                   const ExhaustiveConstVectorView& f_grid,
                                   const AtmPoint& atm,
                                   const zeeman::pol pol) {
  const Numeric H         = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  const Numeric dH_dmag_v = atm.mag[1] / H;

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    dz[i]            = -shp.lines[pos[i].line].inv_gd * dH_dmag_v *
            line.z.Splitting(line.qn.val, pol, pos[i].iz);
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dH(dcut, dz);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dH(dcut, dz, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dH(dz, f);
        });
  }
}

//! Sets dshape and dz
void ComputeData::dmag_w_core_calc(const band_shape& shp,
                                   const band_data& bnd,
                                   const ExhaustiveConstVectorView& f_grid,
                                   const AtmPoint& atm,
                                   const zeeman::pol pol) {
  const Numeric H         = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  const Numeric dH_dmag_w = atm.mag[2] / H;

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    dz[i]            = -shp.lines[pos[i].line].inv_gd * dH_dmag_w *
            line.z.Splitting(line.qn.val, pol, pos[i].iz);
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dH(dcut, dz);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dH(dcut, dz, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dH(dz, f);
        });
  }
}

void compute_derivative(PropmatVectorView dpm,
                        StokvecVectorView dsv,
                        ComputeData& com_data,
                        const ExhaustiveConstVectorView& f_grid,
                        const QuantumIdentifier& qid,
                        const band_shape& shape,
                        const band_data& bnd,
                        const AtmPoint& atm,
                        const zeeman::pol pol,
                        const AtmKey& key) {
  using enum AtmKey;
  switch (key) {
    case t:
      com_data.dt_core_calc(qid, shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dscl[i] * com_data.shape[i].first +
                                    com_data.scl[i] * com_data.dshape[i].first);
        const auto dsv_pm =
            zeeman::scale(com_data.npm,
                          com_data.dscl[i] * com_data.shape[i].second +
                              com_data.scl[i] * com_data.dshape[i].second);
        dsv[i] += {dsv_pm.A(), dsv_pm.B(), dsv_pm.C(), dsv_pm.D()};
      }
      break;
    case p:
      ARTS_USER_ERROR("Not implemented, pressure derivative");
      break;
    case mag_u:
      com_data.dmag_u_core_calc(shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_du,
                                com_data.scl[i] * com_data.shape[i].first,
                                com_data.scl[i] * com_data.dshape[i].first);
        const auto dsv_pm =
            zeeman::scale(com_data.npm,
                          com_data.dnpm_du,
                          com_data.scl[i] * com_data.shape[i].second,
                          com_data.scl[i] * com_data.dshape[i].second);
        dsv[i] += {dsv_pm.A(), dsv_pm.B(), dsv_pm.C(), dsv_pm.D()};
      }
      break;
    case mag_v:
      com_data.dmag_v_core_calc(shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_dv,
                                com_data.scl[i] * com_data.shape[i].first,
                                com_data.scl[i] * com_data.dshape[i].first);
        const auto dsv_pm =
            zeeman::scale(com_data.npm,
                          com_data.dnpm_dv,
                          com_data.scl[i] * com_data.shape[i].second,
                          com_data.scl[i] * com_data.dshape[i].second);
        dsv[i] += {dsv_pm.A(), dsv_pm.B(), dsv_pm.C(), dsv_pm.D()};
      }
      break;
    case mag_w:
      com_data.dmag_w_core_calc(shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_dw,
                                com_data.scl[i] * com_data.shape[i].first,
                                com_data.scl[i] * com_data.dshape[i].first);
        const auto dsv_pm =
            zeeman::scale(com_data.npm,
                          com_data.dnpm_dw,
                          com_data.scl[i] * com_data.shape[i].second,
                          com_data.scl[i] * com_data.dshape[i].second);
        dsv[i] += {dsv_pm.A(), dsv_pm.B(), dsv_pm.C(), dsv_pm.D()};
      }
      break;
    case wind_u:
    case wind_v:
    case wind_w:
      com_data.df_core_calc(shape, bnd, f_grid, atm);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dscl[i] * com_data.shape[i].first +
                                    com_data.scl[i] * com_data.dshape[i].first);
        const auto dsv_pm =
            zeeman::scale(com_data.npm,
                          com_data.dscl[i] * com_data.shape[i].second +
                              com_data.scl[i] * com_data.dshape[i].second);
        dsv[i] += {dsv_pm.A(), dsv_pm.B(), dsv_pm.C(), dsv_pm.D()};
      }
      break;
  }
}

void compute_derivative(PropmatVectorView,
                        StokvecVectorView,
                        ComputeData&,
                        const ExhaustiveConstVectorView&,
                        const QuantumIdentifier& qid,
                        const band_shape&,
                        const band_data&,
                        const AtmPoint&,
                        const zeeman::pol,
                        const SpeciesIsotope& deriv_spec) {
  ARTS_USER_ERROR_IF(deriv_spec == qid.Isotopologue(), "Not supported")
}

void compute_derivative(PropmatVectorView,
                        StokvecVectorView,
                        ComputeData&,
                        const ExhaustiveConstVectorView&,
                        const QuantumIdentifier&,
                        const band_shape&,
                        const band_data&,
                        const AtmPoint&,
                        const zeeman::pol,
                        const SpeciesEnum&) {
  ARTS_USER_ERROR("Not supported")
}

void compute_derivative(PropmatVectorView,
                        StokvecVectorView,
                        ComputeData&,
                        const ExhaustiveConstVectorView&,
                        const QuantumIdentifier&,
                        const band_shape&,
                        const band_data&,
                        const AtmPoint&,
                        const zeeman::pol,
                        const line_key&) {
  ARTS_USER_ERROR("Not supported")
}

void compute_derivative(PropmatVectorView,
                        StokvecVectorView,
                        ComputeData&,
                        const ExhaustiveConstVectorView&,
                        const SpeciesIsotope&,
                        const band_shape&,
                        const band_data&,
                        const AtmPoint&,
                        const zeeman::pol,
                        const auto&) {}

std::ostream& operator<<(std::ostream& os, const ComputeData& cd) {
  for (auto line : cd.lines) {
    os << line.f0 << ' ' << line.k << ' ' << line.e_ratio << ' ' << line.inv_gd
       << ' ' << line.z_imag << '\n';
  }
  return os;
}

void calculate(PropmatVectorView pm,
               StokvecVectorView sv,
               matpack::matpack_view<Propmat, 2, false, true> dpm,
               matpack::matpack_view<Stokvec, 2, false, true> dsv,
               ComputeData& com_data,
               const ExhaustiveConstVectorView& f_grid,
               const JacobianTargets& jacobian_targets,
               const QuantumIdentifier& bnd_qid,
               const band_data& bnd,
               const AtmPoint& atm,
               const zeeman::pol pol,
               const bool no_negative_absorption) {
  ARTS_USER_ERROR_IF(bnd.size() != 1, "Only for single lines per ID")

  if (std::ranges::all_of(com_data.npm, [](auto& n) { return n == 0; })) return;

  const Index nf = f_grid.size();
  if (nf == 0) return;

  const SpeciesIsotope spec = bnd_qid.Isotopologue();
  const Numeric fmin        = f_grid.front();
  const Numeric fmax        = f_grid.back();

  ARTS_ASSERT(jacobian_targets.target_count() ==
                  static_cast<Size>(dpm.nrows()) and
              nf == dpm.ncols())
  ARTS_ASSERT(nf == pm.nelem())

  band_shape_helper(
      com_data.lines, com_data.pos, bnd_qid, bnd, atm, fmin, fmax, pol);

  if (com_data.lines.empty()) return;

  //! Not const to save lines for reuse
  band_shape shape{std::move(com_data.lines), bnd.get_cutoff_frequency()};

  com_data.core_calc(shape, bnd, f_grid);

  for (Index i = 0; i < nf; i++) {
    const auto F = com_data.scl[i] * com_data.shape[i].first;
    const auto N = com_data.scl[i] * com_data.shape[i].second;
    if (no_negative_absorption and F.real() < 0) continue;

    pm[i] += zeeman::scale(com_data.npm, F);

    const auto srcvec  = zeeman::scale(com_data.npm, N);
    sv[i]             += {srcvec.A(), srcvec.B(), srcvec.C(), srcvec.D()};
  }

  for (auto& atm_target : jacobian_targets.atm()) {
    std::visit(
        [&](auto& target) {
          compute_derivative(dpm.as_slice(atm_target.target_pos),
                             dsv.as_slice(atm_target.target_pos),
                             com_data,
                             f_grid,
                             spec,
                             shape,
                             bnd,
                             atm,
                             pol,
                             target);
        },
        atm_target.type);
  }

  for (auto& line_target : jacobian_targets.line()) {
    if (line_target.type.band == bnd_qid) {
      compute_derivative(dpm.as_slice(line_target.target_pos),
                         dsv.as_slice(line_target.target_pos),
                         com_data,
                         f_grid,
                         spec,
                         shape,
                         bnd,
                         atm,
                         pol,
                         line_target.type);
    }
  }

  com_data.lines = std::move(shape.lines);
}
}  // namespace lbl::voigt::nlte
