#include "lbl_lineshape_voigt.h"

#include <new_jacobian.h>
#include <partfun.h>
#include <physics_funcs.h>
#include <sorting.h>

#include <cmath>
#include <cstdio>
#include <limits>

#include <Faddeeva/Faddeeva.hh>

namespace lbl::voigt::lte {
Complex line_strength_calc(const Numeric inv_gd,
                           const SpeciesIsotopeRecord& spec,
                           const line& line,
                           const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const Complex lm{1 + line.ls.G(atm), -line.ls.Y(atm)};

  return Constant::inv_sqrt_pi * inv_gd * atm[spec] * atm[spec.spec] * lm * s;
}

Complex dline_strength_calc_dVMR(const Numeric inv_gd,
                                 const Numeric f0,
                                 const SpeciesIsotopeRecord& spec,
                                 const Species::Species target_spec,
                                 const line& line,
                                 const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const Complex lm{1 + line.ls.G(atm), -line.ls.Y(atm)};

  const Numeric dG = line.ls.dG_dVMR(atm, target_spec);
  const Numeric dY = line.ls.dY_dVMR(atm, target_spec);
  const Numeric dD0 = line.ls.dD0_dVMR(atm, target_spec);
  const Numeric dDV = line.ls.dDV_dVMR(atm, target_spec);
  const Numeric df0 = dD0 + dDV;
  const Complex dlm = {dG, -dY};
  const Numeric r = atm[spec];
  const Numeric x = atm[target_spec];

  if (target_spec == spec.spec) {
    return -Constant::inv_sqrt_pi * inv_gd * r * s *
           (x * (df0 / f0) * lm - (x * dlm + lm));
  }
  return -Constant::inv_sqrt_pi * inv_gd * r * s * x * ((df0 / f0) * lm - dlm);
}

Complex dline_strength_calc_dT(const SpeciesIsotopeRecord& spec,
                               const line& line,
                               const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));
  const auto ds = line.ds_dT(atm.temperature,
                             PartitionFunctions::Q(atm.temperature, spec),
                             PartitionFunctions::dQdT(atm.temperature, spec));

  const Complex lm{1 + line.ls.G(atm), -line.ls.Y(atm)};
  const Complex dlm{line.ls.dG_dT(atm), -line.ls.dY_dT(atm)};

  const auto n = atm[spec] * atm[spec.spec];

  return n * (dlm * s + lm * ds);
}

Complex line_strength_calc(const Numeric inv_gd,
                           const SpeciesIsotopeRecord& spec,
                           const line& line,
                           const AtmPoint& atm,
                           const Size ispec) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const auto& ls = line.ls.single_models[ispec];

  const Complex lm{1 + ls.G(line.ls.T0, atm.temperature, atm.pressure),
                   -ls.Y(line.ls.T0, atm.temperature, atm.pressure)};

  const auto n = atm[spec] * atm[spec.spec] * atm[ls.species];

  return Constant::inv_sqrt_pi * inv_gd * n * lm * s;
}

Complex dline_strength_calc_dVMR(const SpeciesIsotopeRecord& spec,
                                 const Species::Species target_spec,
                                 const line& line,
                                 const AtmPoint& atm,
                                 const Size ispec) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const auto& ls = line.ls.single_models[ispec];

  const Complex lm{1 + ls.G(line.ls.T0, atm.temperature, atm.pressure),
                   -ls.Y(line.ls.T0, atm.temperature, atm.pressure)};

  const auto dn =
      (spec.spec == target_spec ? atm[spec.spec] * atm[ls.species] : 0.0) +
      (spec.spec == ls.species ? atm[spec] * atm[spec.spec] : 0.0);

  return dn * lm * s;
}

Complex dline_strength_calc_dT(const SpeciesIsotopeRecord& spec,
                               const line& line,
                               const AtmPoint& atm,
                               const Size ispec) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));
  const auto ds = line.ds_dT(atm.temperature,
                             PartitionFunctions::Q(atm.temperature, spec),
                             PartitionFunctions::dQdT(atm.temperature, spec));

  const auto& ls = line.ls.single_models[ispec];

  const Complex lm{1 + ls.G(line.ls.T0, atm.temperature, atm.pressure),
                   -ls.Y(line.ls.T0, atm.temperature, atm.pressure)};
  const Complex dlm{ls.dG_dT(line.ls.T0, atm.temperature, atm.pressure),
                    -ls.dY_dT(line.ls.T0, atm.temperature, atm.pressure)};

  const auto n = atm[spec] * atm[spec.spec] * atm[ls.species];

  return n * (dlm * s + lm * ds);
}

Numeric line_center_calc(const line& line, const AtmPoint& atm) {
  return line.f0 + line.ls.D0(atm) + line.ls.DV(atm);
}

Numeric dline_center_calc_dT(const line& line, const AtmPoint& atm) {
  return line.ls.dD0_dT(atm) + line.ls.dDV_dT(atm);
}

Numeric dline_center_calc_dVMR(const line& line,
                               const Species::Species spec,
                               const AtmPoint& atm) {
  return line.ls.dD0_dVMR(atm, spec) + line.ls.dDV_dVMR(atm, spec);
}

Numeric line_center_calc(const line& line, const AtmPoint& atm, Size ispec) {
  const auto& ls = line.ls.single_models[ispec];
  return line.f0 + ls.D0(line.ls.T0, atm.temperature, atm.pressure) +
         ls.DV(line.ls.T0, atm.temperature, atm.pressure);
}

Numeric dline_center_calc_dT(const line& line,
                             const AtmPoint& atm,
                             Size ispec) {
  const auto& ls = line.ls.single_models[ispec];
  return ls.dD0_dT(line.ls.T0, atm.temperature, atm.pressure) +
         ls.dDV_dT(line.ls.T0, atm.temperature, atm.pressure);
}

Numeric dline_center_calc_dVMR(const line& line,
                               const Species::Species spec,
                               const AtmPoint& atm,
                               Size ispec) {
  const auto& ls = line.ls.single_models[ispec];
  return ls.species == spec
             ? ls.D0(line.ls.T0, atm.temperature, atm.pressure) +
                   ls.DV(line.ls.T0, atm.temperature, atm.pressure)
             : 0;
}

Numeric scaled_gd(const Numeric T, const Numeric mass, const Numeric f0) {
  constexpr auto c = Constant::doppler_broadening_const_squared;
  return std::sqrt(c * T / mass) * f0;
}

single_shape::single_shape(const SpeciesIsotopeRecord& spec,
                           const line& line,
                           const AtmPoint& atm)
    : f0(line_center_calc(line, atm)),
      inv_gd(1.0 / scaled_gd(atm.temperature, spec.mass, f0)),
      z_imag(line.ls.G0(atm) * inv_gd),
      s(line_strength_calc(inv_gd, spec, line, atm)) {}

single_shape::single_shape(const SpeciesIsotopeRecord& spec,
                           const line& line,
                           const AtmPoint& atm,
                           const Size ispec)
    : f0(line_center_calc(line, atm, ispec)),
      inv_gd(1.0 / scaled_gd(atm.temperature, spec.mass, f0)),
      z_imag(line.ls.single_models[ispec].G0(
                 line.ls.T0, atm.temperature, atm.pressure) *
             inv_gd),
      s(line_strength_calc(inv_gd, spec, line, atm, ispec)) {}

void single_shape::as_zeeman(const line& line,
                             const Numeric H,
                             zeeman::pol pol,
                             Index iz) {
  s *= line.z.Strength(line.qn.val, pol, iz);
  f0 += H * line.z.Splitting(line.qn.val, pol, iz);
}

Complex single_shape::F(const Complex z_) noexcept { return Faddeeva::w(z_); }

Complex single_shape::F(const Numeric f) const noexcept { return F(z(f)); }

Complex single_shape::operator()(const Numeric f) const noexcept {
  return s * F(f);
}

Complex single_shape::dF(const Numeric f) const noexcept {
  const Complex z_ = z(f);
  return dF(z_, F(z_));
}

Complex single_shape::dF(const Complex z_, const Complex F_) noexcept {
  /*! FIXME: We should use a proper algorithm here.  This produces
   *         no errors in tests, but the actual derivative form is
   *         analytically known.  Its numerical instability, however,
   *         makes it completely useless.
   *
   *         The analytical form is:
   *
   *           dF = -2 * z * F(z) + 2 * i / sqrt(pi)
  */
  const Complex dz = 1e-6 * z_;
  const Complex F_2 = Faddeeva::w(z_ + dz);
  return (F_2 - F_) / dz;
}

single_shape::zFdF::zFdF(const Complex z_) noexcept
    : z{z_}, F{single_shape::F(z_)}, dF{single_shape::dF(z_, F)} {}

single_shape::zFdF single_shape::all(const Numeric f) const noexcept {
  return z(f);
}

Complex single_shape::df(const Numeric f) const noexcept { return s * dF(f); }

Complex single_shape::df0(const Complex ds_df0,
                          const Complex dz_df0,
                          const Numeric f) const noexcept {
  const auto [z_, F_, dF_] = all(f);
  return ds_df0 * F_ + s * dz_df0 * dF_;
}

Complex single_shape::dDV(const Complex dz_dDV,
                          const Numeric f) const noexcept {
  return s * dz_dDV * dF(f);
}

Complex single_shape::dD0(const Complex dz_dD0,
                          const Numeric f) const noexcept {
  return s * dz_dD0 * dF(f);
}

Complex single_shape::dG0(const Complex dz_dG0,
                          const Numeric f) const noexcept {
  return s * dz_dG0 * dF(f);
}

Complex single_shape::dH(const Complex dz_dH, const Numeric f) const noexcept {
  return s * dz_dH * dF(f);
}

Complex single_shape::dVMR(const Complex ds_dVMR,
                           const Complex dz_dVMR,
                           const Numeric dz_dVMR_fac,
                           const Numeric f) const noexcept {
  const auto [z_, F_, dF_] = all(f);
  return ds_dVMR * F_ + s * (dz_dVMR + dz_dVMR_fac * z_) * dF_;
}

Complex single_shape::dT(const Complex ds_dT,
                         const Complex dz_dT,
                         const Numeric f) const noexcept {
  const auto [z_, F_, dF_] = all(f);
  return ds_dT * F_ + dz_dT * dF_;
  //FIXME: invGD factor of F_ is missing
}

Complex single_shape::da(const Complex ds_da, const Numeric f) const noexcept {
  return ds_da * F(f);
}

Complex single_shape::de0(const Complex ds_de0,
                          const Numeric f) const noexcept {
  return ds_de0 * F(f);
}

Complex single_shape::dG(const Complex ds_dG, const Numeric f) const noexcept {
  return ds_dG * F(f);
}

Complex single_shape::dY(const Complex ds_dY, const Numeric f) const noexcept {
  return ds_dY * F(f);
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
                      single_shape&& s,
                      const line& line,
                      const Numeric H,
                      const Size spec,
                      const zeeman::pol pol) {
  const auto line_nr = static_cast<Size>(pos.size() ? pos.back().line + 1 : 0);

  if (pol == zeeman::pol::no) {
    lines.emplace_back(std::forward<single_shape>(s));
    pos.emplace_back(line_nr, spec);
  } else {
    const auto nz = line.z.size(line.qn.val, pol);
    for (Index iz = 0; iz < nz; iz++) {
      lines.push_back(s);
      lines.back().as_zeeman(line, H, pol, iz);
      pos.emplace_back(line_nr, spec, static_cast<Size>(iz));
    }
  }
}

void lines_push_back(std::vector<single_shape>& lines,
                     std::vector<line_pos>& pos,
                     const SpeciesIsotopeRecord& spec,
                     const line& line,
                     const AtmPoint& atm,
                     const zeeman::pol pol) {
  const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  if (line.ls.one_by_one) {
    for (Size i = 0; i < line.ls.single_models.size(); ++i) {
      zeeman_push_back(
          lines, pos, single_shape{spec, line, atm, i}, line, H, i, pol);
    }
  } else {
    zeeman_push_back(lines,
                     pos,
                     single_shape{spec, line, atm},
                     line,
                     H,
                     std::numeric_limits<Size>::max(),
                     pol);
  }
}

void band_shape_helper(std::vector<single_shape>& lines,
                       std::vector<line_pos>& pos,
                       const SpeciesIsotopeRecord& spec,
                       const band_data& bnd,
                       const AtmPoint& atm,
                       const Numeric fmin,
                       const Numeric fmax,
                       const zeeman::pol pol) {
  lines.resize(0);
  pos.resize(0);
  lines.reserve(count_lines(bnd, pol));
  pos.reserve(lines.capacity());

  using enum CutoffType;
  switch (bnd.cutoff) {
    case None:
      for (auto& line : bnd) {
        lines_push_back(lines, pos, spec, line, atm, pol);
      }
      break;
    case ByLine: {
      const auto active_lines = bnd.active_lines(fmin, fmax);
      for (auto& line : active_lines) {
        lines_push_back(lines, pos, spec, line, atm, pol);
      }
    } break;
    case FINAL:
      ARTS_USER_ERROR("Bad state")
  }

  bubble_sort_by(
      [&](const Size l1, const Size l2) { return lines[l1].f0 < lines[l2].f0; },
      lines,
      pos);
}

band_shape::band_shape(std::vector<single_shape>&& ls,
                       const Numeric cut) noexcept
    : lines(std::move(ls)), cutoff(cut) {}

Complex band_shape::operator()(const Numeric f) const noexcept {
  return std::transform_reduce(
      lines.begin(), lines.end(), Complex{}, std::plus<>{}, [f](auto& ls) {
        return ls(f);
      });
}

Complex band_shape::df(const Numeric f) const {
  return std::transform_reduce(
      lines.begin(), lines.end(), Complex{}, std::plus<>{}, [f](auto& ls) {
        return ls.df(f);
      });
}

Complex band_shape::dH(const ExhaustiveConstComplexVectorView& dz_dH,
                       const Numeric f) const {
  ARTS_ASSERT(static_cast<Size>(dz_dH.size()) == lines.size())

  return std::transform_reduce(lines.begin(),
                               lines.end(),
                               dz_dH.begin(),
                               Complex{},
                               std::plus<>{},
                               [f](auto& ls, auto& d) { return ls.dH(d, f); });
}

Complex band_shape::dT(const ExhaustiveConstComplexVectorView& ds_dT,
                       const ExhaustiveConstComplexVectorView& dz_dT,
                       const Numeric f) const {
  ARTS_ASSERT(ds_dT.size() == dz_dT.size())
  ARTS_ASSERT(static_cast<Size>(ds_dT.size()) == lines.size())

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i = 0; i < lines.size(); ++i) {
    out += lines[i].dT(ds_dT[i], dz_dT[i], f);
  }

  return out;
}

Complex band_shape::dVMR(const ExhaustiveConstComplexVectorView& ds_dVMR,
                         const ExhaustiveConstComplexVectorView& dz_dVMR,
                         const ExhaustiveConstVectorView& dz_dVMR_fac,
                         const Numeric f) const {
  ARTS_ASSERT(ds_dVMR.size() == dz_dVMR.size())
  ARTS_ASSERT(static_cast<Size>(ds_dVMR.size()) == lines.size())

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i = 0; i < lines.size(); ++i) {
    out += lines[i].dVMR(ds_dVMR[i], dz_dVMR[i], dz_dVMR_fac[i], f);
  }

  return out;
}

Complex band_shape::df0(const ExhaustiveConstComplexVectorView ds_df0,
                        const ExhaustiveConstComplexVectorView dz_df0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].df0(ds_df0[i], dz_df0[i], f);
  }

  return out;
}

Complex band_shape::da(const ExhaustiveConstComplexVectorView ds_da,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].da(ds_da[i], f);
  }

  return out;
}

Complex band_shape::de0(const ExhaustiveConstComplexVectorView ds_de0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].de0(ds_de0[i], f);
  }

  return out;
}

Complex band_shape::dDV(const ExhaustiveConstComplexVectorView dz_dDV,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dDV(dz_dDV[i], f);
  }

  return out;
}

Complex band_shape::dD0(const ExhaustiveConstComplexVectorView dz_dD0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dD0(dz_dD0[i], f);
  }

  return out;
}

Complex band_shape::dG0(const ExhaustiveConstComplexVectorView dz_dG0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dG0(dz_dG0[i], f);
  }

  return out;
}

Complex band_shape::dY(const ExhaustiveConstComplexVectorView ds_dY,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dY(ds_dY[i], f);
  }

  return out;
}

Complex band_shape::dG(const ExhaustiveConstComplexVectorView ds_dG,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dG(ds_dG[i], f);
  }

  return out;
}

Complex band_shape::operator()(const ExhaustiveConstComplexVectorView& cut,
                               const Numeric f) const {
  const auto [s, cs] = frequency_spans(cutoff, f, lines, cut);
  return std::transform_reduce(s.begin(),
                               s.end(),
                               cs.begin(),
                               Complex{},
                               std::plus<>{},
                               [f](auto& ls, auto& c) { return ls(f) - c; });
}

void band_shape::operator()(ExhaustiveComplexVectorView cut) const {
  std::transform(
      lines.begin(),
      lines.end(),
      cut.begin(),
      [cutoff_freq = cutoff](auto& ls) { return ls(ls.f0 + cutoff_freq); });
}

Complex band_shape::df(const ExhaustiveConstComplexVectorView& cut,
                       const Numeric f) const {
  const auto [s, cs] = frequency_spans(cutoff, f, lines, cut);
  return std::transform_reduce(s.begin(),
                               s.end(),
                               cs.begin(),
                               Complex{},
                               std::plus<>{},
                               [f](auto& ls, auto& c) { return ls.df(f) - c; });
}

void band_shape::df(ExhaustiveComplexVectorView cut) const {
  std::transform(
      lines.begin(),
      lines.end(),
      cut.begin(),
      [cutoff_freq = cutoff](auto& ls) { return ls.df(ls.f0 + cutoff_freq); });
}

Complex band_shape::dH(const ExhaustiveConstComplexVectorView& cut,
                       const ExhaustiveConstComplexVectorView& dz_dH,
                       const Numeric f) const {
  ARTS_ASSERT(static_cast<Size>(dz_dH.size()) == lines.size())

  const auto [s, cs, dH] = frequency_spans(cutoff, f, lines, cut, dz_dH);

  Complex out{};  //! Fixme, use zip in C++ 23...
  for (Size i = 0; i < s.size(); ++i) {
    out += s[i].dH(dH[i], f) - cs[i];
  }

  return out;
}

void band_shape::dH(ExhaustiveComplexVectorView cut,
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

Complex band_shape::dT(const ExhaustiveConstComplexVectorView& cut,
                       const ExhaustiveConstComplexVectorView& ds_dT,
                       const ExhaustiveConstComplexVectorView& dz_dT,
                       const Numeric f) const {
  ARTS_ASSERT(ds_dT.size() == dz_dT.size())
  ARTS_ASSERT(static_cast<Size>(ds_dT.size()) == lines.size())

  Complex out{};  //! Fixme, use zip in C++ 23...

  const auto [s, cs, ds, dz] =
      frequency_spans(cutoff, f, lines, cut, ds_dT, dz_dT);

  for (Size i = 0; i < s.size(); ++i) {
    out += s[i].dT(ds[i], dz[i], f) - cs[i];
  }

  return out;
}

void band_shape::dT(ExhaustiveComplexVectorView cut,
                    const ExhaustiveConstComplexVectorView& ds_dT,
                    const ExhaustiveConstComplexVectorView& dz_dT) const {
  ARTS_ASSERT(ds_dT.size() == dz_dT.size())
  ARTS_ASSERT(static_cast<Size>(ds_dT.size()) == lines.size())

  for (Size i = 0; i < lines.size(); ++i) {
    cut[i] += lines[i].dT(ds_dT[i], dz_dT[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dVMR(const ExhaustiveConstComplexVectorView& cut,
                         const ExhaustiveConstComplexVectorView& ds_dVMR,
                         const ExhaustiveConstComplexVectorView& dz_dVMR,
                         const ExhaustiveConstVectorView& dz_dVMR_fac,
                         const Numeric f) const {
  ARTS_ASSERT(ds_dVMR.size() == dz_dVMR.size())
  ARTS_ASSERT(static_cast<Size>(ds_dVMR.size()) == lines.size())

  Complex out{};  //! Fixme, use zip in C++ 23...

  const auto [s, cs, ds, dz, dzf] =
      frequency_spans(cutoff, f, lines, cut, ds_dVMR, dz_dVMR, dz_dVMR_fac);

  for (Size i = 0; i < s.size(); ++i) {
    out += s[i].dVMR(ds[i], dz[i], dzf[i], f) - cs[i];
  }

  return out;
}

void band_shape::dVMR(ExhaustiveComplexVectorView cut,
                      const ExhaustiveConstComplexVectorView& ds_dVMR,
                      const ExhaustiveConstComplexVectorView& dz_dVMR,
                      const ExhaustiveConstVectorView& dz_dVMR_fac) const {
  ARTS_ASSERT(ds_dVMR.size() == dz_dVMR.size())
  ARTS_ASSERT(static_cast<Size>(ds_dVMR.size()) == lines.size())

  for (Size i = 0; i < lines.size(); ++i) {
    cut[i] += lines[i].dVMR(
        ds_dVMR[i], dz_dVMR[i], dz_dVMR_fac[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::df0(const ExhaustiveConstComplexVectorView& cut,
                        const ExhaustiveConstComplexVectorView ds_df0,
                        const ExhaustiveConstComplexVectorView dz_df0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  const auto [s, cs, ds, dz] =
      frequency_spans(cutoff, f, lines, cut, ds_df0, dz_df0);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].df0(ds[i], dz[i], f) - cs[i];
  }

  return out;
}

void band_shape::df0(ExhaustiveComplexVectorView cut,
                     const ExhaustiveConstComplexVectorView ds_df0,
                     const ExhaustiveConstComplexVectorView dz_df0,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].df0(ds_df0[i], dz_df0[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::da(const ExhaustiveConstComplexVectorView& cut,
                       const ExhaustiveConstComplexVectorView ds_da,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_da);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].da(ds[i], f) - cs[i];
  }

  return out;
}

void band_shape::da(ExhaustiveComplexVectorView cut,
                    const ExhaustiveComplexVectorView ds_da,
                    const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].da(ds_da[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::de0(const ExhaustiveConstComplexVectorView& cut,
                        const ExhaustiveConstComplexVectorView ds_de0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_de0);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].de0(ds[i], f) - cs[i];
  }

  return out;
}

void band_shape::de0(ExhaustiveComplexVectorView cut,
                     const ExhaustiveComplexVectorView ds_de0,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].de0(ds_de0[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dDV(const ExhaustiveConstComplexVectorView& cut,
                        const ExhaustiveConstComplexVectorView dz_dDV,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  const auto [s, cs, dz] = frequency_spans(cutoff, f, lines, cut, dz_dDV);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dDV(dz[i], f) - cs[i];
  }

  return out;
}

void band_shape::dDV(ExhaustiveComplexVectorView cut,
                     const ExhaustiveComplexVectorView dz_dDV,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].dDV(dz_dDV[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dD0(const ExhaustiveConstComplexVectorView& cut,
                        const ExhaustiveConstComplexVectorView dz_dD0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  const auto [s, cs, dz] = frequency_spans(cutoff, f, lines, cut, dz_dD0);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dD0(dz[i], f) - cs[i];
  }

  return out;
}

void band_shape::dD0(ExhaustiveComplexVectorView cut,
                     const ExhaustiveComplexVectorView dz_dD0,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].dD0(dz_dD0[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dG0(const ExhaustiveConstComplexVectorView& cut,
                        const ExhaustiveConstComplexVectorView dz_dG0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  const auto [s, cs, dz] = frequency_spans(cutoff, f, lines, cut, dz_dG0);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dG0(dz[i], f) - cs[i];
  }

  return out;
}

void band_shape::dG0(ExhaustiveComplexVectorView cut,
                     const ExhaustiveComplexVectorView dz_dG0,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].dG0(dz_dG0[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dY(const ExhaustiveConstComplexVectorView& cut,
                       const ExhaustiveConstComplexVectorView ds_dY,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_dY);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dY(ds[i], f) - cs[i];
  }

  return out;
}

void band_shape::dY(ExhaustiveComplexVectorView cut,
                    const ExhaustiveComplexVectorView ds_dY,
                    const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].dY(ds_dY[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dG(const ExhaustiveConstComplexVectorView& cut,
                       const ExhaustiveConstComplexVectorView ds_dG,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_dG);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dG(ds[i], f) - cs[i];
  }

  return out;
}

void band_shape::dG(ExhaustiveComplexVectorView cut,
                    const ExhaustiveComplexVectorView ds_dG,
                    const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].dG(ds_dG[i], lines[i].f0 + cutoff);
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
  using Constant::h, Constant::k;

  std::transform(
      f_grid.begin(),
      f_grid.end(),
      scl.begin(),
      [N = number_density(atm.pressure, atm.temperature), T = atm.temperature](
          auto f) { return -N * f * std::expm1(-h * f / (k * T)); });

  update_zeeman(los, atm.mag, pol);
}

//! Sizes cut, dcut, dz, ds; sets shape
void ComputeData::core_calc(const band_shape& shp,
                            const band_data& bnd,
                            const ExhaustiveConstVectorView& f_grid) {
  cut.resize(shp.size());
  dz.resize(shp.size());
  dz_fac.resize(shp.size());
  ds.resize(shp.size());
  dcut.resize(shp.size());
  filter.reserve(shp.size());

  if (bnd.cutoff != CutoffType::None) {
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
void ComputeData::dt_core_calc(const SpeciesIsotopeRecord& spec,
                               const band_shape& shp,
                               const band_data& bnd,
                               const ExhaustiveConstVectorView& f_grid,
                               const AtmPoint& atm,
                               const zeeman::pol pol) {
  using Constant::h, Constant::k;

  std::transform(f_grid.begin(),
                 f_grid.end(),
                 dscl.begin(),
                 [N = number_density(atm.pressure, atm.temperature),
                  dN = dnumber_density_dt(atm.pressure, atm.temperature),
                  T = atm.temperature](auto f) {
                   return -f *
                          (N * f * h * exp(-h * f / (k * T)) / (T * T * k) +
                           dN * std::expm1(-h * f / (k * T)));
                 });

  const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  for (Size i = 0; i < pos.size(); i++) {
    const bool single_line = pos[i].spec == std::numeric_limits<Size>::max();
    const Numeric& inv_gd = shp.lines[i].inv_gd;
    const auto& line = bnd.lines[pos[i].line];
    const Numeric dinv_gd_dT_ratio = -0.5 / atm.temperature;
    ds[i] = dinv_gd_dT_ratio * shp.lines[i].s +
            Constant::inv_sqrt_pi * inv_gd *
                line.z.Strength(line.qn.val, pol, pos[i].iz) *
                (single_line
                     ? dline_strength_calc_dT(spec, line, atm)
                     : dline_strength_calc_dT(spec, line, atm, pos[i].spec));
    dz[i] =
        dinv_gd_dT_ratio *
            Complex{-inv_gd * shp.lines[i].f0, shp.lines[i].z_imag} +
        inv_gd *
            Complex{
                -(single_line ? dline_center_calc_dT(line, atm)
                              : dline_center_calc_dT(line, atm, pos[i].spec)) -
                    H * line.z.Splitting(line.qn.val, pol, pos[i].iz),
                (single_line ? line.ls.dG0_dT(atm)
                             : line.ls.single_models[pos[i].spec].dG0_dT(
                                   line.ls.T0, atm.temperature, atm.pressure))};
  }

  if (bnd.cutoff != CutoffType::None) {
    shp.dT(dcut, ds, dz);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dT(dcut, ds, dz, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dT(ds, dz, f);
        });
  }
}

//! Sets dshape and dscl
void ComputeData::df_core_calc(const band_shape& shp,
                               const band_data& bnd,
                               const ExhaustiveConstVectorView& f_grid,
                               const AtmPoint& atm) {
  using Constant::h, Constant::k;

  std::transform(f_grid.begin(),
                 f_grid.end(),
                 dscl.begin(),
                 [N = number_density(atm.pressure, atm.temperature),
                  T = atm.temperature](auto f) {
                   return N * (f * h * std::exp(-f * h / (T * k)) / (T * k) -
                               std::expm1(-f * h / (T * k)));
                 });

  if (bnd.cutoff != CutoffType::None) {
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
  const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  const Numeric dH_dmag_u = atm.mag[0] / H;

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    dz[i] = -shp.lines[pos[i].line].inv_gd * dH_dmag_u *
            line.z.Splitting(line.qn.val, pol, pos[i].iz);
  }

  if (bnd.cutoff != CutoffType::None) {
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
  const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  const Numeric dH_dmag_v = atm.mag[1] / H;

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    dz[i] = -shp.lines[pos[i].line].inv_gd * dH_dmag_v *
            line.z.Splitting(line.qn.val, pol, pos[i].iz);
  }

  if (bnd.cutoff != CutoffType::None) {
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
  const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  const Numeric dH_dmag_w = atm.mag[2] / H;

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    dz[i] = -shp.lines[pos[i].line].inv_gd * dH_dmag_w *
            line.z.Splitting(line.qn.val, pol, pos[i].iz);
  }

  if (bnd.cutoff != CutoffType::None) {
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

//! Sets ds and dz and dcut and dshape
void ComputeData::dVMR_core_calc(const SpeciesIsotopeRecord& spec,
                                 const band_shape& shp,
                                 const band_data& bnd,
                                 const ExhaustiveConstVectorView& f_grid,
                                 const AtmPoint& atm,
                                 const zeeman::pol pol,
                                 const Species::Species target_spec) {
  using Constant::h, Constant::k;

  for (Size i = 0; i < pos.size(); i++) {
    const Numeric& inv_gd = shp.lines[i].inv_gd;
    const Numeric& f0 = shp.lines[i].f0;
    const auto& line = bnd.lines[pos[i].line];

    if (pos[i].spec == std::numeric_limits<Size>::max()) {
      dz_fac[i] = -(line.ls.dD0_dVMR(atm, target_spec) +
                    line.ls.dDV_dVMR(atm, target_spec)) /
                  f0;

      ds[i] =
          line.z.Strength(line.qn.val, pol, pos[i].iz) *
          dline_strength_calc_dVMR(inv_gd, f0, spec, target_spec, line, atm);

      dz[i] = inv_gd * Complex{-dline_center_calc_dVMR(line, target_spec, atm),
                               line.ls.dG0_dVMR(atm, target_spec)};
    } else {
      ARTS_ASSERT(false, "Not implemented")
    }
  }

  if (bnd.cutoff != CutoffType::None) {
    shp.dVMR(dcut, ds, dz, dz_fac);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dVMR(dcut, ds, dz, dz_fac, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dVMR(ds, dz, dz_fac, f);
        });
  }
}

void ComputeData::set_filter(const line_key& key, bool check_spec) {
  if (key.line == filtered_line and
      (not check_spec or (check_spec and key.spec == filtered_spec)))
    return;

  filtered_spec = check_spec ? key.spec : std::numeric_limits<Size>::max();
  filtered_line = key.line;
  filter.resize(0);

  if (not check_spec) {
    for (Size i = 0; i < pos.size(); i++) {
      if (pos[i].line == key.line) filter.push_back(i);
    }
  } else {
    for (Size i = 0; i < pos.size(); i++) {
      if (pos[i].line == key.line and pos[i].spec == key.spec)
        filter.push_back(i);
    }
  }
}

//! Sets dshape and ds and dz and dcut and dshape
void ComputeData::df0_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ExhaustiveConstVectorView& f_grid,
                                const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key, false);
  const Numeric ds_df0_ratio =
      bnd.lines[pos[filter.front()].line].ds_df0_s_ratio();

  for (Size i : filter) {
    const auto& line = shp.lines[i];
    const Numeric& inv_gd = line.inv_gd;
    const Numeric dinv_gd_df0_ratio = -1 / line.f0;
    ds[i] = dinv_gd_df0_ratio * line.s + ds_df0_ratio * line.s;
    dz[i] =
        dinv_gd_df0_ratio * Complex{-inv_gd * line.f0, line.z_imag} - inv_gd;
  }

  if (bnd.cutoff != CutoffType::None) {
    shp.df0(dcut, ds, dz, filter);
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.df0(dcut, ds, dz, f_grid[i], filter);
    }
  } else {
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.df0(ds, dz, f_grid[i], filter);
    }
  }
}

//! Sets dshape and ds and dcut and dshape
void ComputeData::de0_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ExhaustiveConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key, false);

  for (Size i : filter) {
    const Numeric ds_de0_ratio =
        bnd.lines[pos[i].line].ds_de0_s_ratio(atm.temperature);
    ds[i] = ds_de0_ratio * shp.lines[i].s;
  }

  if (bnd.cutoff != CutoffType::None) {
    shp.de0(dcut, ds, filter);
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.de0(dcut, ds, f_grid[i], filter);
    }
  } else {
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.de0(ds, f_grid[i], filter);
    }
  }
}

//! Sets dshape and ds and dcut and dshape
void ComputeData::da_core_calc(const band_shape& shp,
                               const band_data& bnd,
                               const ExhaustiveConstVectorView& f_grid,
                               const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key, false);

  for (Size i : filter) {
    const Numeric ds_da_ratio = 1.0 / bnd.lines[pos[i].line].a;
    ds[i] = ds_da_ratio * shp.lines[i].s;
  }

  if (bnd.cutoff != CutoffType::None) {
    shp.da(dcut, ds, filter);
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.da(dcut, ds, f_grid[i], filter);
    }
  } else {
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.da(ds, f_grid[i], filter);
    }
  }
}

//! Sets dshape and dz and dcut and dshape
void ComputeData::dG0_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ExhaustiveConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key, true);

  using enum temperature::coefficient;
  switch (key.ls_coeff) {
    case X0:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        dz[i] = Complex(
            0,
            shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG0_dX0(
                                      ls.T0, atm.temperature, atm.pressure));
      }
      break;
    case X1:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        dz[i] = Complex(
            0,
            shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG0_dX1(
                                      ls.T0, atm.temperature, atm.pressure));
      }
      break;
    case X2:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        dz[i] = Complex(
            0,
            shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG0_dX2(
                                      ls.T0, atm.temperature, atm.pressure));
      }
      break;
    case X3:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        dz[i] = Complex(
            0,
            shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG0_dX3(
                                      ls.T0, atm.temperature, atm.pressure));
      }
      break;
    case FINAL:;
  }

  if (bnd.cutoff != CutoffType::None) {
    shp.dG0(dcut, dz, filter);
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dG0(dcut, dz, f_grid[i], filter);
    }
  } else {
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dG0(dz, f_grid[i], filter);
    }
  }
}

//! Sets dshape and dz and dcut and dshape
void ComputeData::dD0_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ExhaustiveConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key, true);

  using enum temperature::coefficient;
  switch (key.ls_coeff) {
    case X0:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        dz[i] = -lines[i].inv_gd * ls.single_models[pos[i].spec].dD0_dX0(
                                       ls.T0, atm.temperature, atm.pressure);
      }
      break;
    case X1:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        dz[i] =
            -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dD0_dX1(
                                       ls.T0, atm.temperature, atm.pressure);
      }
      break;
    case X2:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        dz[i] =
            -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dD0_dX2(
                                       ls.T0, atm.temperature, atm.pressure);
      }
      break;
    case X3:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        dz[i] =
            -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dD0_dX3(
                                       ls.T0, atm.temperature, atm.pressure);
      }
      break;
    case FINAL:;
  }

  if (bnd.cutoff != CutoffType::None) {
    shp.dD0(dcut, dz, filter);
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dD0(dcut, dz, f_grid[i], filter);
    }
  } else {
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dD0(dz, f_grid[i], filter);
    }
  }
}

//! Sets dshape and ds and dcut and dshape
void ComputeData::dY_core_calc(const band_shape& shp,
                               const band_data& bnd,
                               const ExhaustiveConstVectorView& f_grid,
                               const AtmPoint& atm,
                               const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key, true);

  using enum temperature::coefficient;
  switch (key.ls_coeff) {
    case X0:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        ds[i] = Complex{
            0,
            -lines[i].inv_gd * ls.single_models[pos[i].spec].dY_dX0(
                                   ls.T0, atm.temperature, atm.pressure)};
      }
      break;
    case X1:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        ds[i] = Complex{
            0,
            -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dY_dX1(
                                       ls.T0, atm.temperature, atm.pressure)};
      }
      break;
    case X2:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        ds[i] = Complex{
            0,
            -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dY_dX2(
                                       ls.T0, atm.temperature, atm.pressure)};
      }
      break;
    case X3:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        ds[i] = Complex{
            0,
            -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dY_dX3(
                                       ls.T0, atm.temperature, atm.pressure)};
      }
      break;
    case FINAL:;
  }

  if (bnd.cutoff != CutoffType::None) {
    shp.dY(dcut, dz, filter);
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dY(dcut, ds, f_grid[i], filter);
    }
  } else {
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dY(ds, f_grid[i], filter);
    }
  }
}

//! Sets dshape and ds and dcut and dshape
void ComputeData::dG_core_calc(const band_shape& shp,
                               const band_data& bnd,
                               const ExhaustiveConstVectorView& f_grid,
                               const AtmPoint& atm,
                               const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key, true);

  using enum temperature::coefficient;
  switch (key.ls_coeff) {
    case X0:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        ds[i] = lines[i].inv_gd * ls.single_models[pos[i].spec].dG_dX0(
                                      ls.T0, atm.temperature, atm.pressure);
      }
      break;
    case X1:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        ds[i] = shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG_dX1(
                                          ls.T0, atm.temperature, atm.pressure);
      }
      break;
    case X2:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        ds[i] = shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG_dX2(
                                          ls.T0, atm.temperature, atm.pressure);
      }
      break;
    case X3:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        ds[i] = shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG_dX3(
                                          ls.T0, atm.temperature, atm.pressure);
      }
      break;
    case FINAL:;
  }

  if (bnd.cutoff != CutoffType::None) {
    shp.dG(dcut, dz, filter);
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dG(dcut, ds, f_grid[i], filter);
    }
  } else {
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dG(ds, f_grid[i], filter);
    }
  }
}

//! Sets dshape and dz and dcut and dshape
void ComputeData::dDV_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ExhaustiveConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key, true);

  using enum temperature::coefficient;
  switch (key.ls_coeff) {
    case X0:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        dz[i] = -lines[i].inv_gd * ls.single_models[pos[i].spec].dDV_dX0(
                                       ls.T0, atm.temperature, atm.pressure);
      }
      break;
    case X1:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        dz[i] =
            -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dDV_dX1(
                                       ls.T0, atm.temperature, atm.pressure);
      }
      break;
    case X2:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        dz[i] =
            -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dDV_dX2(
                                       ls.T0, atm.temperature, atm.pressure);
      }
      break;
    case X3:
      for (Size i : filter) {
        const auto& ls = bnd.lines[pos[i].line].ls;
        dz[i] =
            -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dDV_dX3(
                                       ls.T0, atm.temperature, atm.pressure);
      }
      break;
    case FINAL:;
  }

  if (bnd.cutoff != CutoffType::None) {
    shp.dDV(dcut, dz, filter);
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dDV(dcut, dz, f_grid[i], filter);
    }
  } else {
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dDV(dz, f_grid[i], filter);
    }
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const ExhaustiveConstVectorView& f_grid,
                        const SpeciesIsotopeRecord& spec,
                        const band_shape& shape,
                        const band_data& bnd,
                        const AtmPoint& atm,
                        const zeeman::pol pol,
                        const Atm::Key& key) {
  using enum Atm::Key;
  switch (key) {
    case t:
      com_data.dt_core_calc(spec, shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dscl[i] * com_data.shape[i] +
                                    com_data.scl[i] * com_data.dshape[i]);
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
                                com_data.scl[i] * com_data.shape[i],
                                com_data.scl[i] * com_data.dshape[i]);
      }
      break;
    case mag_v:
      com_data.dmag_v_core_calc(shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_dv,
                                com_data.scl[i] * com_data.shape[i],
                                com_data.scl[i] * com_data.dshape[i]);
      }
      break;
    case mag_w:
      com_data.dmag_w_core_calc(shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_dw,
                                com_data.scl[i] * com_data.shape[i],
                                com_data.scl[i] * com_data.dshape[i]);
      }
      break;
    case wind_u:
    case wind_v:
    case wind_w:
      com_data.df_core_calc(shape, bnd, f_grid, atm);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dscl[i] * com_data.shape[i] +
                                    com_data.scl[i] * com_data.dshape[i]);
      }
      break;
    case FINAL:
      break;
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const ExhaustiveConstVectorView& f_grid,
                        const SpeciesIsotopeRecord& spec,
                        const band_shape&,
                        const band_data&,
                        const AtmPoint& atm,
                        const zeeman::pol,
                        const SpeciesIsotopeRecord& deriv_spec) {
  if (deriv_spec != spec) return;

  const Numeric isorat = atm[spec];

  ARTS_USER_ERROR_IF(
      isorat == 0,
      "Does not support 0 for isotopologue ratios (may be added upon request)")

  for (Index i = 0; i < f_grid.size(); i++) {
    dpm[i] += zeeman::scale(com_data.npm,
                            com_data.scl[i] * com_data.shape[i] / isorat);
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const ExhaustiveConstVectorView& f_grid,
                        const SpeciesIsotopeRecord& spec,
                        const band_shape& shape,
                        const band_data& bnd,
                        const AtmPoint& atm,
                        const zeeman::pol pol,
                        const Species::Species& deriv_spec) {
  com_data.dVMR_core_calc(spec, shape, bnd, f_grid, atm, pol, deriv_spec);
  for (Index i = 0; i < f_grid.size(); i++) {
    dpm[i] += zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const ExhaustiveConstVectorView& f_grid,
                        const band_shape& shape,
                        const band_data& bnd,
                        const AtmPoint& atm,
                        const zeeman::pol,
                        const line_key& deriv) {
  switch (deriv.var) {
    case variable::f0:
      com_data.df0_core_calc(shape, bnd, f_grid, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case variable::e0:
      com_data.de0_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case variable::a:
      com_data.da_core_calc(shape, bnd, f_grid, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case variable::FINAL:
      break;
  }

  switch (deriv.ls_var) {
    case line_shape::variable::G0:
      com_data.dG0_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case line_shape::variable::D0:
      com_data.dD0_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case line_shape::variable::G2:
      return;
    case line_shape::variable::D2:
      return;
    case line_shape::variable::FVC:
      return;
    case line_shape::variable::ETA:
      return;
    case line_shape::variable::Y:
      com_data.dY_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case line_shape::variable::G:
      com_data.dG_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case line_shape::variable::DV:
      com_data.dDV_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case line_shape::variable::FINAL:
      break;
  }
}

void compute_derivative(PropmatVectorView,
                        ComputeData&,
                        const ExhaustiveConstVectorView&,
                        const SpeciesIsotopeRecord&,
                        const band_shape&,
                        const band_data&,
                        const AtmPoint&,
                        const zeeman::pol,
                        const auto&) {}

void calculate(PropmatVectorView pm,
               matpack::matpack_view<Propmat, 2, false, true> dpm,
               ComputeData& com_data,
               const ExhaustiveConstVectorView& f_grid,
               const JacobianTargets& jacobian_targets,
               const QuantumIdentifier& bnd_qid,
               const band_data& bnd,
               const AtmPoint& atm,
               const zeeman::pol pol) {
  const Index nf = f_grid.size();
  if (nf == 0) return;

  const SpeciesIsotopeRecord spec = bnd_qid.Isotopologue();
  const Numeric fmin = f_grid.front();
  const Numeric fmax = f_grid.back();

  ARTS_ASSERT(jacobian_targets.target_count() ==
                  static_cast<Size>(dpm.nrows()) and
              nf == dpm.ncols())
  ARTS_ASSERT(nf == pm.nelem())

  //! Reuse lines and positions if possible
  band_shape_helper(
      com_data.lines, com_data.pos, spec, bnd, atm, fmin, fmax, pol);

  if (com_data.lines.empty()) return;

  //! Not const to save lines for reuse
  band_shape shape{std::move(com_data.lines), bnd.get_cutoff_frequency()};

  com_data.core_calc(shape, bnd, f_grid);

  for (Index i = 0; i < nf; i++) {
    pm[i] += zeeman::scale(com_data.npm, com_data.scl[i] * com_data.shape[i]);
  }

  for (auto& atm_target : jacobian_targets.atm()) {
    std::visit(
        [&](auto& target) {
          compute_derivative(dpm.as_slice(atm_target.target_pos),
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
                         com_data,
                         f_grid,
                         shape,
                         bnd,
                         atm,
                         pol,
                         line_target.type);
    }
  }

  com_data.lines = std::move(shape.lines);
}
}  // namespace lbl::voigt::lte
